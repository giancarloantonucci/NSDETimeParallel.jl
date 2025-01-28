"MPI implementation of Parareal."
function parareal_mpi!(
    cache::PararealCache, solution::PararealSolution, problem::AbstractInitialValueProblem, parareal::Parareal;
    directory::String="results", saveiterates::Bool=false)
    
    # Extract components
    @↓ finesolver, coarsesolver = parareal
    @↓ N, K = parareal.parameters
    @↓ weights, ψ, ϵ = parareal.tolerance
    @↓ makeGs, U, T, F, G = cache
    @↓ errors = solution

    # MPI set-up
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    root = 0

    # Ensure directory exists (only rank 0 creates it)
    if rank == root && !isdir(directory)
        mkpath(directory)
    end
    # Synchronize so all ranks see the directory
    MPI.Barrier(comm)

    if rank == root
        # Initialization
        F[1] = G[1] = U[1]

        # Coarse run (serial)
        for n = 1:N-1
            chunkproblem = copy(problem, U[n], T[n], T[n+1])
            chunkcoarsesolution = coarsesolver(chunkproblem)
            G[n+1] = chunkcoarsesolution(T[n+1])
            if !makeGs[n+1]
                U[n+1] = copy(G[n+1])
            end
        end
    end

    # Broadcast initial values to all processes
    T = MPI.bcast(T, root, comm)

    # Main loop
    isconverged = false
    chunkfinesolution = nothing # needed for savelastiterate
    for k = 1:K

        # TODO: Don't need to bcast the whole of U
        U = MPI.bcast(U, root, comm)
        
        # Save U and T at current iteration
        # if rank == root && saveiterates
        #     iterates[k].U .= U
        #     iterates[k].T .= T
        # end

        # Fine run (parallelised with MPI): Rank 0 receives fine results from ranks k..N-1
        n = rank
        if k ≤ n ≤ N
            chunkproblem = copy(problem, U[n], T[n], T[n+1])
            chunkfinesolution = finesolver(chunkproblem)
            if n < N
                U_n_plus_1 = chunkfinesolution(T[n+1])
                MPI.send(U_n_plus_1, comm, dest=root, tag=n)
            end
            # Save chunks at current iteration
            if saveiterates
                filename = joinpath(directory, "iter_$(k)_chunk_$(n).jls")
                open(filename, "w") do file
                    local_data = (chunk_n = chunkfinesolution,)
                    serialize(file, local_data)
                end
            end
        end
        if rank == root
            for n = k:N-1
                F[n+1] = MPI.recv(comm, source=n, tag=n)
            end
        end

        if rank == root
            # Update weights of error function
            update!(weights, U, F)

            # Check convergence
            errors[k] = ψ(cache, k, weights)
            if errors[k] ≤ ϵ
                resize!(errors, k) # self-updates inside solution
                isconverged = true
            end
        end

        isconverged = MPI.bcast(isconverged, root, comm)
        if isconverged
            if rank != root
                n = rank
                filename = joinpath(directory, "lastiter_chunk_$(n).jls")
                open(filename, "w") do file
                    local_data = (chunk_n = chunkfinesolution,)
                    serialize(file, local_data)
                end
            end
            break
        end

        if rank == root
            # Correction step (serial)
            for n = k:N-1
                chunkproblem = copy(problem, U[n], T[n], T[n+1])
                chunkcoarsesolution = coarsesolver(chunkproblem)
                v = chunkcoarsesolution(T[n+1])
                U[n+1] = v + F[n+1] - G[n+1]
                G[n+1] = v
            end
        end
    end

    MPI.Barrier(comm)
    
    if rank == root
        collect!(solution; directory)
    end

    return solution
end
