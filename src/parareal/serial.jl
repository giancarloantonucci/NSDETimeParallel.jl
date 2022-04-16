"Serial implementation of Parareal."
function parareal_serial!(cache::PararealCache, solution::PararealSolution, problem::AbstractInitialValueProblem, parareal::Parareal)
    @↓ U, F, G, T = cache
    @↓ errors, alliterates, saveiterates = solution
    @↓ finesolver, coarsolver, P, K = parareal
    @↓ weights, ψ, ϵ = parareal.tolerance
    # coarse run (serial)
    # ...outside
    # main loop
    F[1] = U[1]
    for k in 1:K
        # fine run (parallelisable)
        for n = k:P
            chunkproblem = subproblemof(problem, U[n], T[n], T[n+1])
            chunksolution = finesolver(chunkproblem)
            solution[n] = chunksolution
            n < P ? F[n+1] = chunksolution.u[end] : nothing
        end
        # save iterates
        if saveiterates
            for n = 1:k-1
                alliterates[k][n] = alliterates[k-1][n]
            end
            for n = k:P
                alliterates[k][n] = solution[n]
            end
        end
        # check convergence
        update!(weights, U, F)
        errors[k] = ψ(cache, solution, k, weights)
        if errors[k] ≤ ϵ
            resize!(errors, k)
            if saveiterates
                resize!(alliterates, k)
            end
            break
        end
        # serial update
        for n = k:P-1
            chunkproblem = subproblemof(problem, U[n], T[n], T[n+1])
            chunksolution = coarsolver(chunkproblem)
            U[n+1] = chunksolution.u[end] + F[n+1] - G[n+1]
            G[n+1] = chunksolution.u[end]
        end
    end
    return solution
end
