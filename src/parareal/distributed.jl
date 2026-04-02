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

    F[1] = G[1] = U[1]
    
    for n = 1:N-1
        chunkproblem = copy(problem, U[n], T[n], T[n+1])
        chunkcoarsesolution = coarsesolver(chunkproblem)
        G[n+1] = chunkcoarsesolution(T[n+1])
        if makeGs[n+1]
            U[n+1] = copy(G[n+1])
        end
    end

    for k = 1:K
        # Pure mathematical mapping without global variables or disk I/O
        function finesolve(n)
            local_problem = copy(problem, U[n], T[n], T[n+1])
            local_solution = finesolver(local_problem)
            return local_solution(T[n+1]), local_solution
        end
        
        # Dispatch to workers
        fineresults = pmap(finesolve, WorkerPool(workers()[k:N]), k:N)
        
        # CRITICAL FIX: Loop up to N to ensure all chunks are saved.
        # F[n+1] is only updated if n < N.
        for n = k:N
            final_state, local_sol = fineresults[n-k+1]
            solution[n] = local_sol
            
            if n < N
                F[n+1] = final_state
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
            chunkcoarsesolution = coarsesolver(chunkproblem)
            v = chunkcoarsesolution(T[n+1])
            
            U[n+1] = v + F[n+1] - G[n+1]
            G[n+1] = v
        end
    end

    # Serialise the collected solutions from the master process memory
    if !nocollect
        for n = 1:N
            if isassigned(solution.lastiterate.chunks, n)
                filename = joinpath(directory, "lastiter_chunk_$(n).jls")
                open(filename, "w") do file
                    serialize(file, (chunk_n = solution[n],))
                end
            end
        end
    end

    return solution
end
