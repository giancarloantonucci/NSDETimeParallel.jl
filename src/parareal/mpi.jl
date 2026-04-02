# src/parareal/mpi.jl

"MPI implementation of Parareal."
function parareal_mpi!(
    cache::PararealCache, solution::PararealSolution,
    problem::AbstractInitialValueProblem, parareal::Parareal;
    directory::String, saveiterates::Bool, nocollect::Bool)
    
    @↓ finesolver, coarsesolver = parareal
    @↓ N, K = parareal.parameters
    @↓ weights, ψ, ϵ = parareal.tolerance
    @↓ makeGs, U, T, F, G = cache
    @↓ errors = solution

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    size = MPI.Comm_size(comm)
    root = 0

    # CRITICAL FIX: Align 0-indexed MPI ranks with 1-indexed Julia chunks
    n = rank + 1 

    if rank == root && !isdir(directory)
        mkpath(directory)
    end
    MPI.Barrier(comm)

    if rank == root
        F[1] = G[1] = U[1]
        # Coarse run (strictly serial on master)
        for m = 1:N-1
            chunkproblem = copy(problem, U[m], T[m], T[m+1])
            chunkcoarsesolution = coarsesolver(chunkproblem)
            G[m+1] = chunkcoarsesolution(T[m+1])
            if makeGs[m+1]
                U[m+1] = copy(G[m+1])
            end
        end
    end

    T = MPI.bcast(T, root, comm)

    isconverged = false
    local_fine_chunk = nothing

    for k = 1:K
        # 1. Distribute specific starting points
        if rank == root
            for dest_rank = 1:size-1
                w = dest_rank + 1 # chunk index for the worker
                if w >= k
                    MPI.Send(U[w], comm; dest=dest_rank, tag=w)
                end
            end
        elseif n >= k
            if !isassigned(U, n) U[n] = similar(U[1]) end
            MPI.Recv!(U[n], comm; source=root, tag=n)
        end

        # 2. Fine run (parallel execution across all ranks, including root)
        if n >= k && n <= N
            chunkproblem = copy(problem, U[n], T[n], T[n+1])
            local_fine_chunk = finesolver(chunkproblem)
            
            if saveiterates
                filename = joinpath(directory, "iter_$(k)_chunk_$(n).jls")
                open(filename, "w") do file
                    serialize(file, (chunk_n = local_fine_chunk,))
                end
            end
            
            if n < N
                U_n_plus_1 = local_fine_chunk(T[n+1])
                # Worker sends result back to root
                if rank != root
                    MPI.Send(U_n_plus_1, comm; dest=root, tag=n+N)
                else
                    F[n+1] = U_n_plus_1 # Root handles its own chunk locally
                end
            end
        end

        # 3. Master collects results and calculates error
        if rank == root
            for src_rank = 1:size-1
                w = src_rank + 1
                if w >= k && w < N
                    if !isassigned(F, w+1) F[w+1] = similar(U[1]) end
                    MPI.Recv!(F[w+1], comm; source=src_rank, tag=w+N)
                end
            end
            
            update!(weights, U, F)
            errors[k] = ψ(cache, k, weights)
            
            if errors[k] ≤ ϵ
                resize!(errors, k)
                isconverged = true
            end
        end

        # Broadcast convergence status to all ranks
        isconverged = MPI.bcast(isconverged, root, comm)
        if isconverged
            break
        end

        # 4. Correction step (serial on master)
        if rank == root
            for m = k:N-1
                chunkproblem = copy(problem, U[m], T[m], T[m+1])
                chunkcoarsesolution = coarsesolver(chunkproblem)
                v = chunkcoarsesolution(T[m+1])
                
                U[m+1] = v + F[m+1] - G[m+1]
                G[m+1] = v
            end
        end
    end

    # Serialise cleanly at the end
    if local_fine_chunk !== nothing
        filename = joinpath(directory, "lastiter_chunk_$(n).jls")
        open(filename, "w") do file
            serialize(file, (chunk_n = local_fine_chunk,))
        end
    end

    MPI.Barrier(comm)
    
    if rank == root && !nocollect
        collect!(solution; directory)
    end

    return solution
end
