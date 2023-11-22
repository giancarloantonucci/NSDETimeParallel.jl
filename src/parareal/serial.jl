"Serial implementation of Parareal."
function parareal_serial!(cache::PararealCache, solution::PararealSolution, problem::AbstractInitialValueProblem, parareal::Parareal)
    @↓ skips, F, G = cache
    @↓ errors, iterates = solution
    @↓ U, T = solution.lastiterate
    @↓ finesolver, coarsesolver, saveiterates = parareal
    @↓ N, K = parareal.parameters
    @↓ weights, ψ, ϵ = parareal.tolerance

    # coarse run (serial)
    for n = 1:N
        if skips[n] # `skips[1] == true` always
            G[n] = U[n]
        else
            chunkproblem = copy(problem, U[n-1], T[n-1], T[n])
            chunkcoarsesolution = coarsesolver(chunkproblem)
            G[n] = chunkcoarsesolution(T[n])
        end
    end

    # initialization
    F[1] = U[1]

    # main loop
    for k in 1:K

        # fine run (parallelisable)
        for n = k:N
            chunkproblem = copy(problem, U[n], T[n], T[n+1])
            chunksolution = finesolver(chunkproblem)
            solution[n] = chunksolution
            if n < N
                F[n+1] = chunksolution(T[n+1])
            end
        end

        # save iterates
        if saveiterates
            iterates[k].U .= U
            iterates[k].T .= T
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
            chunkproblem = copy(problem, U[n], T[n], T[n+1])
            chunksolution = coarsesolver(chunkproblem)
            v = chunksolution(T[n+1])
            U[n+1] = v + F[n+1] - G[n+1]
            G[n+1] = v
        end
        
    end

    return solution
end
