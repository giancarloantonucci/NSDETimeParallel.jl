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

    # Align 0-indexed MPI ranks with 1-indexed Julia chunks
    n = rank + 1 

    if rank == root && !isdir(directory)
        mkpath(directory)
    end
    MPI.Barrier(comm)

    # 1. Pre-allocate zero-allocation caches
    dummy_chunk = copy(problem, U[1], T[1], T[2])
    
    # Every rank handles exactly one time chunk per iteration, so it only needs ONE fine cache.
    local_fine_cache = NSDEBase.initialize_cache(dummy_chunk, finesolver)
    local_fine_sol   = NSDEBase.initialize_solution(dummy_chunk, finesolver)
    
    # Only the master coordinates the coarse solves serially, so it only needs ONE coarse cache.
    if rank == root
        master_coarse_cache = NSDEBase.initialize_cache(dummy_chunk, coarsesolver)
        master_coarse_sol   = NSDEBase.initialize_solution(dummy_chunk, coarsesolver)
    end

    # Allocate the inner buffers for U, G, F to prevent UndefRefErrors
    for m = 2:N
        if !isassigned(U, m) U[m] = similar(U[1]) end
        if !isassigned(G, m) G[m] = similar(U[1]) end
        if !isassigned(F, m) F[m] = similar(U[1]) end
    end

    if rank == root
        F[1] = G[1] = U[1]
        # Coarse run (strictly serial on master)
        for m = 1:N-1
            chunkproblem = copy(problem, U[m], T[m], T[m+1])
            NSDEBase.solve!(master_coarse_cache, master_coarse_sol, chunkproblem, coarsesolver)
            
            # Must copy to prevent alias bleeding when master_coarse_sol is reused
            copyto!(G[m+1], master_coarse_sol(T[m+1]))
            if makeGs[m+1]
                copyto!(U[m+1], G[m+1])
            end
        end
    end

    T = MPI.bcast(T, root, comm)

    isconverged = false

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
            MPI.Recv!(U[n], comm; source=root, tag=n)
        end

        # 2. Fine run (parallel execution across all ranks, including root)
        if n >= k && n <= N
            chunkproblem = copy(problem, U[n], T[n], T[n+1])
            NSDEBase.solve!(local_fine_cache, local_fine_sol, chunkproblem, finesolver)
            
            if saveiterates
                filename = joinpath(directory, "iter_$(k)_chunk_$(n).jls")
                open(filename, "w") do file
                    # Must deepcopy because local_fine_sol will be overwritten
                    serialize(file, (chunk_n = deepcopy(local_fine_sol),))
                end
            end
            
            if n < N
                U_n_plus_1 = local_fine_sol(T[n+1])
                # Worker sends result back to root
                if rank != root
                    MPI.Send(U_n_plus_1, comm; dest=root, tag=n+N)
                else
                    copyto!(F[n+1], U_n_plus_1) # Root handles its own chunk locally
                end
            end
        end

        # 3. Master collects results and calculates error
        if rank == root
            for src_rank = 1:size-1
                w = src_rank + 1
                if w >= k && w < N
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
                NSDEBase.solve!(master_coarse_cache, master_coarse_sol, chunkproblem, coarsesolver)
                
                v = master_coarse_sol(T[m+1])
                
                # Rebind U to protect U_previous in the error function ψ₁
                U[m+1] = v + F[m+1] - G[m+1]
                
                # Explicit copyto! for G to prevent alias bleeding from reused coarse sol
                copyto!(G[m+1], v)
            end
        end
    end

    # Serialise cleanly at the end
    if n <= N
        filename = joinpath(directory, "lastiter_chunk_$(n).jls")
        open(filename, "w") do file
            serialize(file, (chunk_n = deepcopy(local_fine_sol),))
        end
    end

    MPI.Barrier(comm)
    
    if rank == root && !nocollect
        collect!(solution; directory)
    end

    return solution
end
