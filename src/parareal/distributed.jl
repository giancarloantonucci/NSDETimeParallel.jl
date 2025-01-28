"Distributed implementation of Parareal."
function parareal_distributed!(cache::PararealCache, solution::PararealSolution, problem::AbstractInitialValueProblem, parareal::Parareal)
    @↓ skips, F, G = cache
    @↓ errors, iterates = solution
    @↓ U, T = solution.lastiterate
    @↓ u0, (t0, tN) ← tspan = problem
    @↓ finesolver, coarsesolver, saveiterates = parareal
    @↓ N, K = parareal.parameters
    @↓ weights, ψ, ϵ = parareal.tolerance

    # Initialization
    for n = 1:N+1
        T[n] = (N - n + 1) / N * t0 + (n - 1) / N * tN # stable sum
    end
    F[1] = G[1] = U[1] = u0

    # Coarse run (serial)
    for n = 1:N-1
        chunkproblem = copy(problem, U[n], T[n], T[n+1])
        chunkcoarsesolution = coarsesolver(chunkproblem)
        G[n+1] = chunkcoarsesolution(T[n+1])
        if !skips[n+1]
            U[n+1] = G[n+1]
        end
    end

    # Main loop
    for k in 1:K

        # Fine run (parallelised with pmap)
        # docs: @distributed for many simple operations / pmap for a few complex operations
        function finesolve(worker_rank)
            n = worker_rank
            chunkproblem = copy(problem, U[n], T[n], T[n+1])
            global chunkfinesolution # save full chunk solution locally for later retrieval
            chunkfinesolution = finesolver(chunkproblem)
            Uₚ = chunkfinesolution(T[n+1])
            return Uₚ
        end
        fineresults = pmap(finesolve, WorkerPool(workers()[k:N]), k:N)

        for n = k:N-1
            F[n+1] = fineresults[n-k+1]
        end

        # Save iterates
        if saveiterates
            iterates[k].U .= U
            iterates[k].T .= T
            for n = 1:k-1
                iterates[k][n] = iterates[k-1][n]
            end
            for n = k:N
                iterates[k][n] = @fetchfrom workers()[n] NSDETimeParallel.chunkfinesolution
            end
        end

        # Update weights of error function
        update!(weights, U, F)

        # Check convergence
        errors[k] = ψ(cache, solution, k, weights)
        if errors[k] ≤ ϵ
            resize!(errors, k) # self-updates inside solution
            if saveiterates
                resize!(iterates, k) # self-updates inside solution
            end
            break
        end

        # Correction step (serial)
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
