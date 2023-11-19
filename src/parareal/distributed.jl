"Distributed implementation of Parareal."
function parareal_distributed!(cache::PararealCache, solution::PararealSolution, problem::AbstractInitialValueProblem, parareal::Parareal)
    @↓ skips, U, F, G, T = cache
    @↓ errors, iterates = solution
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

        # fine run (parallelised with pmap)
        function finesolve(worker_rank)
            n = worker_rank
            chunkproblem = copy(problem, U[n], T[n], T[n+1])
            global chunkfinesolution # save full chunk solution for later retrieval
            chunkfinesolution = finesolver(chunkproblem)
            Uₚ = chunkfinesolution(T[n+1])
            return Uₚ
        end
        fineresults = pmap(finesolve, WorkerPool(workers()[k:N]), k:N)

        for n = k:N-1
            F[n+1] = fineresults[n-k+1]
        end

        # save iterates
        if saveiterates
            for n = 1:k-1
                iterates[k][n] = iterates[k-1][n]
            end
            for n = k:N
                iterates[k][n] = @fetchfrom workers()[n] NSDETimeParallel.chunkfinesolution
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
            chunkcoarsesolution = coarsesolver(chunkproblem)
            v = chunkcoarsesolution(T[n+1])
            U[n+1] = v + F[n+1] - G[n+1]
            G[n+1] = v
        end

    end

    return solution
end
