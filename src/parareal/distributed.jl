# src/parareal/distributed.jl

"Distributed implementation of Parareal."
function parareal_distributed!(
    cache::PararealCache, solution::PararealSolution,
    problem::AbstractInitialValueProblem, parareal::Parareal;
    directory::String, saveiterates::Bool, nocollect::Bool)
    
    @↓ finesolver, coarsesolver = parareal
    @↓ N, K = parareal.parameters
    @↓ weights, ψ, ϵ = parareal.tolerance
    @↓ makeGs, U, T, F, G = cache
    @↓ errors = solution

    if !isdir(directory)
        mkpath(directory)
    end

    # 1. Pre-allocate zero-allocation cache for the coarse solver (Master Node only)
    dummy_chunk = copy(problem, U[1], T[1], T[2])
    master_coarse_cache = NSDEBase.initialize_cache(dummy_chunk, coarsesolver)
    master_coarse_sol   = NSDEBase.initialize_solution(dummy_chunk, coarsesolver)

    F[1] = G[1] = U[1]

    # Allocate inner buffers to give copyto! physical memory
    for n = 2:N
        if !isassigned(U, n) U[n] = similar(U[1]) end
        if !isassigned(G, n) G[n] = similar(U[1]) end
        if !isassigned(F, n) F[n] = similar(U[1]) end
    end
    
    # Coarse run (serial on master)
    for n = 1:N-1
        chunkproblem = copy(problem, U[n], T[n], T[n+1])
        NSDEBase.solve!(master_coarse_cache, master_coarse_sol, chunkproblem, coarsesolver)
        
        copyto!(G[n+1], master_coarse_sol(T[n+1]))
        if makeGs[n+1]
            copyto!(U[n+1], G[n+1])
        end
    end

    for k = 1:K
        # Pure mathematical mapping without returning heavy arrays
        function finesolve(n)
            local_problem = copy(problem, U[n], T[n], T[n+1])
            local_solution = finesolver(local_problem)
            
            # Worker writes to disk directly to bypass the IPC bottleneck
            if saveiterates
                filename = joinpath(directory, "iter_$(k)_chunk_$(n).jls")
                open(filename, "w") do file
                    serialize(file, (chunk_n = local_solution,))
                end
            end
            
            filename = joinpath(directory, "lastiter_chunk_$(n).jls")
            open(filename, "w") do file
                serialize(file, (chunk_n = local_solution,))
            end
            
            # ONLY return the lightweight boundary value to the master
            return local_solution(T[n+1])
        end
        
        # Dispatch to workers
        fineresults = pmap(finesolve, WorkerPool(workers()[k:N]), k:N)
        
        # Master updates F with the lightweight boundary values
        for n = k:N
            final_state = fineresults[n-k+1]
            if n < N
                copyto!(F[n+1], final_state)
            end
        end

        update!(weights, U, F)

        errors[k] = ψ(cache, k, weights)
        if errors[k] ≤ ϵ
            resize!(errors, k)
            break
        end

        for n = k:N-1
            chunkproblem = copy(problem, U[n], T[n], T[n+1])
            NSDEBase.solve!(master_coarse_cache, master_coarse_sol, chunkproblem, coarsesolver)
            
            v = master_coarse_sol(T[n+1])
            
            U[n+1] = v + F[n+1] - G[n+1]
            copyto!(G[n+1], v)
        end
    end

    # Master reads the collected solutions from disk back into memory
    if !nocollect
        collect!(solution; directory)
    end

    return solution
end
