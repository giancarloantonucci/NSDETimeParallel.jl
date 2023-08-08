"Distributed implementation of Parareal."
function parareal_distributed!(cache::PararealCache, solution::PararealSolution, problem::AbstractInitialValueProblem, parareal::Parareal)
  @↓ skips, U, F, G, T = cache
  @↓ errors, iterates = solution
  @↓ finesolver, coarsolver, saveiterates = parareal
  @↓ N, K = parareal.parameters
  @↓ weights, ψ, ϵ = parareal.tolerance

  # coarse run (serial)
  for n = 1:N
    if skips[n] # `skips[1] == true` always
      G[n] = U[n]
    else
      chunkproblem = subproblemof(problem, U[n-1], T[n-1], T[n])
      G[n] = coarsolver(chunkproblem)(T[n])
    end
  end

  # initialization
  F[1] = U[1]

  # send fine solver function to procs
  function finesolve(finesolver, chunkproblem)
    chunksolution = finesolver(chunkproblem)
    Uₚ = chunksolution.u[end]
    return Uₚ, chunksolution
  end

  # main loop
  for k in 1:K

    # fine run (parallelised with Distributed.jl)
    tasks = []
    for worker_rank in procs()
      n = worker_rank
      chunkproblem = subproblemof(problem, U[n], T[n], T[n+1])
      push!(tasks, remotecall(finesolve, n, finesolver, chunkproblem))
    end

    fineresults = fetch.(tasks)

    for n = k:N
      solution[n] = fineresults[n][2]
      if n < N
        F[n+1] = fineresults[n][1]
      end
    end

    # save iterates
    if saveiterates
      for n = 1:k-1
        iterates[k][n] = iterates[k-1][n]
      end
      for n = k:N
        iterates[k][n] = solution[n]
      end
    end

    # check convergence
    update!(weights, U, F)
    errors[k] = ψ(cache, solution, k, weights)
    if errors[k] ≤ ϵ
      resize!(errors, k) # self-updates in solution
      if saveiterates
        resize!(iterates, k) # self-updates in solution
      end
      break
    end

    # correction step (serial)
    for n = k:N-1
      chunkproblem = subproblemof(problem, U[n], T[n], T[n+1])
      chunksolution = coarsolver(chunkproblem)
      v = chunksolution(T[n+1])
      U[n+1] = v + F[n+1] - G[n+1]
      G[n+1] = v
    end
  end

  return solution
end
