"Distributed implementation of Parareal."
function parareal_distributed!(
    cache::PararealCache, solution::PararealSolution, problem::AbstractInitialValueProblem, parareal::Parareal;
    directory::String="results", saveiterates::Bool=false)

    # Extract components
    @↓ skips, U, T, F, G = cache
    @↓ errors = solution
    @↓ finesolver, coarsesolver = parareal
    @↓ N, K = parareal.parameters
    @↓ weights, ψ, ϵ = parareal.tolerance

    # Ensure directory exists (only rank 0 creates it)
    if !isdir(directory)
        mkpath(directory)
    end

    # Initialization
    F[1] = G[1] = U[1]
        
    # Coarse run (serial)
    for n = 1:N-1
        chunkproblem = copy(problem, U[n], T[n], T[n+1])
        chunkcoarsesolution = coarsesolver(chunkproblem)
        G[n+1] = chunkcoarsesolution(T[n+1])
        if !skips[n+1]
            U[n+1] = copy(G[n+1])
        end
    end

    # Main loop
    for k = 1:K

        # Fine run (parallelised with pmap)
        # docs: @distributed for many simple operations / pmap for a few complex operations
        function finesolve(worker_rank)
            n = worker_rank
            chunkproblem = copy(problem, U[n], T[n], T[n+1])
            global chunkfinesolution # save full chunk solution locally for later retrieval
            chunkfinesolution = finesolver(chunkproblem)
            # Save chunks at current iteration
            if saveiterates
                filename = joinpath(directory, "iter_$(k)_chunk_$(n).jls")
                open(filename, "w") do file
                    local_data = (chunk_n = chunkfinesolution,)
                    serialize(file, local_data)
                end
            end
            U_n_plus_1 = chunkfinesolution(T[n+1])
            return U_n_plus_1
        end
        fineresults = pmap(finesolve, WorkerPool(workers()[k:N]), k:N)
        for n = k:N-1
            F[n+1] = fineresults[n-k+1]
        end

        # Update weights of error function
        update!(weights, U, F)

        # Check convergence
        errors[k] = ψ(cache, k, weights)
        if errors[k] ≤ ϵ
            resize!(errors, k) # self-updates inside solution
            function saveresults(worker_rank)
                n = worker_rank
                filename = joinpath(directory, "lastiter_chunk_$(n).jls")
                open(filename, "w") do file
                    local_data = (chunk_n = chunkfinesolution,)
                    serialize(file, local_data)
                end
            end
            pmap(saveresults, WorkerPool(workers()[1:N]), 1:N)
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

    collect!(solution; directory)

    return solution
end
