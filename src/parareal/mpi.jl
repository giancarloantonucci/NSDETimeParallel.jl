"MPI implementation of Parareal."
function parareal_mpi!(cache::PararealCache, solution::PararealSolution, problem::AbstractInitialValueProblem, parareal::Parareal)
    
    # MPI set-up
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    size = MPI.Comm_size(comm)
    root = 0
    
    @↓ skips, F, G = cache
    @↓ errors, iterates = solution
    @↓ U, T = solution.lastiterate
    @↓ u0, (t0, tN) ← tspan = problem
    @↓ finesolver, coarsesolver, saveiterates = parareal
    @↓ N, K = parareal.parameters
    @↓ weights, ψ, ϵ = parareal.tolerance

    # Initialization
    if rank == root
        for n = 1:N+1
            T[n] = (N - n + 1) / N * t0 + (n - 1) / N * tN # stable sum
        end
        F[1] = G[1] = U[1] = u0
    end

    # Broadcast initial values to all processes
    T = MPI.bcast(T, root, comm)
    U = MPI.bcast(U, root, comm)

    # Coarse run (serial)
    if rank == root
        for n = 1:N-1
            chunkproblem = copy(problem, U[n], T[n], T[n+1])
            chunkcoarsesolution = coarsesolver(chunkproblem)
            G[n+1] = chunkcoarsesolution(T[n+1])
            if !skips[n+1]
                U[n+1] = G[n+1]
            end
        end
    end
    U = MPI.bcast(U, root, comm)

    # Main loop
    isconverged = false
    for k in 1:K

        # Fine run (parallelised with MPI)
        # function finesolve(worker_rank)
        #     n = worker_rank
        #     chunkproblem = copy(problem, U[n], T[n], T[n+1])
        #     global chunkfinesolution # save full chunk solution locally for later retrieval
        #     chunkfinesolution = finesolver(chunkproblem)
        #     Uₚ = chunkfinesolution(T[n+1])
        #     return Uₚ
        # end
        # fineresults = pmap(finesolve, WorkerPool(workers()[k:N]), k:N)
        local_fine_results = Vector{Any}(undef, 0)
        for n = k:N - 1
            if (n % size == rank)
                chunkproblem = copy(problem, U[n], T[n], T[n+1])
                global chunkfinesolution # save full chunk solution locally for later retrieval
                chunkfinesolution = finesolver(chunkproblem)
                Uₚ = chunkfinesolution(T[n+1])
                push!(local_fine_results, (n, Uₚ))
            end
        end

        # for n = k:N-1
        #     F[n+1] = fineresults[n-k+1]
        # end
        # Gather results at rank 0
        all_fine_results = MPI.gather(local_fine_results, comm, root=root)
        if rank == root
            for results in all_fine_results
                for (n, Uₚ) in results
                    F[n+1] = Uₚ
                end
            end
        end

        # Save iterates
        if saveiterates && rank == root
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
        if rank == root
            update!(weights, U, F)
        end

        # Check convergence
        if rank == root
            errors[k] = ψ(cache, solution, k, weights)
            if errors[k] ≤ ϵ
                resize!(errors, k) # self-updates inside solution
                if saveiterates
                    resize!(iterates, k) # self-updates inside solution
                end
                isconverged = true
            end
        end
        isconverged = MPI.bcast(isconverged, root, comm)
        if isconverged
            break
        end

        # Correction step (serial)
        if rank == root
            for n = k:N-1
                chunkproblem = copy(problem, U[n], T[n], T[n+1])
                chunkcoarsesolution = coarsesolver(chunkproblem)
                v = chunkcoarsesolution(T[n+1])
                U[n+1] = v + F[n+1] - G[n+1]
                G[n+1] = v
            end
        end
        U = MPI.bcast(U, root, comm)
    end

    return solution
end
