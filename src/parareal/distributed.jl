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
        G[n] = coarsesolver(chunkproblem)(T[n])
        end
    end

    # initialization
    F[1] = U[1]

    # send fine solver function to procs (when parallelised with remotecall and fetch)
    # function finesolve(finesolver, chunkproblem, saveiterates)
    #     chunksolution = finesolver(chunkproblem)
    #     Uₚ = chunksolution.u[end]
    #     if saveiterates
    #         return (Uₚ, chunksolution)
    #     else
    #         return Uₚ
    #     end
    # end

    # main loop
    for k in 1:K

        # fine run (parallelised with remotecall and fetch) NB: SLOWER!
        # tasks = Vector{Future}()
        # for worker_rank in k:N
        #     n = worker_rank
        #     chunkproblem = copy(problem, U[n], T[n], T[n+1])
        #     push!(tasks, remotecall(finesolve, n, finesolver, chunkproblem, saveiterates))
        # end
        # fineresults = fetch.(tasks)

        # fine run (parallelised with distributed)
        # combine(destination, source) = append!(destination, source)
        # fineresults = @sync @distributed combine for worker_rank in k:N
        #     n = worker_rank
        #     chunkproblem = copy(problem, U[n], T[n], T[n+1])
        #     chunksolution = finesolver(chunkproblem)
        #     Uₚ = chunksolution(T[n+1])
        #     if saveiterates
        #         return [(Uₚ, chunksolution)]
        #     else
        #         return [Uₚ]
        #     end
        # end

        # fine run (parallelised with pmap)
        function finesolve(worker_rank)
            n = worker_rank
            chunkproblem = copy(problem, U[n], T[n], T[n+1])
            chunksolution = finesolver(chunkproblem)
            Uₚ = chunksolution.u[end]
            if saveiterates
                return (Uₚ, chunksolution)
            else
                return Uₚ
            end
        end
        fineresults = pmap(finesolve, WorkerPool(workers()[k:N]), k:N)

        for n = k:N
            if saveiterates
                solution[n] = fineresults[n-k+1][2]
                if n < N
                    F[n+1] = fineresults[n-k+1][1]
                end
            else
                if n < N
                    F[n+1] = fineresults[n-k+1]
                end
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
            chunkproblem = copy(problem, U[n], T[n], T[n+1])
            chunksolution = coarsesolver(chunkproblem)
            v = chunksolution(T[n+1])
            U[n+1] = v + F[n+1] - G[n+1]
            G[n+1] = v
        end

    end

    return solution
end
