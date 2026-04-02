# src/parareal/threads.jl

"Multi-threaded implementation of Parareal (Shared Memory)."
function parareal_threads!(
    cache::PararealCache, solution::PararealSolution,
    problem::AbstractInitialValueProblem, parareal::Parareal;
    directory::String, saveiterates::Bool, nocollect::Bool)
    
    # Extract components
    @↓ finesolver, coarsesolver = parareal
    @↓ N, K = parareal.parameters
    @↓ weights, ψ, ϵ = parareal.tolerance
    
    # Get U and T from cache
    @↓ makeGs, U, T, F, G = cache 
    @↓ errors = solution

    # 1. Pre-allocate thread-safe caches and solution objects
    dummy_chunk = copy(problem, U[1], T[1], T[2])
    fine_caches = [NSDEBase.initialize_cache(dummy_chunk, finesolver) for _ in 1:N]
    fine_sols   = [NSDEBase.initialize_solution(dummy_chunk, finesolver) for _ in 1:N]

    # Initialization
    F[1] = G[1] = U[1]
    
    # Allocate the inner buffers for U, G, F to prevent UndefRefErrors
    for n = 2:N
        if !isassigned(U, n) U[n] = similar(U[1]) end
        if !isassigned(G, n) G[n] = similar(U[1]) end
        if !isassigned(F, n) F[n] = similar(U[1]) end
    end
    
    # Coarse run (strictly serial)
    for n = 1:N-1
        chunkproblem = copy(problem, U[n], T[n], T[n+1])
        # TODO: For maximum performance, pre-allocate the coarse solver cache too!
        chunkcoarsesolution = coarsesolver(chunkproblem)
        
        copyto!(G[n+1], chunkcoarsesolution(T[n+1])) 
        if makeGs[n+1]
            copyto!(U[n+1], G[n+1])
        end
    end

    # Main loop
    for k = 1:K
        # Fine run (Multi-threaded loop)
        Threads.@threads :dynamic for n = k:N
            chunkproblem = copy(problem, U[n], T[n], T[n+1])
            
            local_cache = fine_caches[n]
            local_sol   = fine_sols[n]
            
            # Zero-allocation solve
            NSDEBase.solve!(local_cache, local_sol, chunkproblem, finesolver)
            
            solution[n] = local_sol 

            if n < N
                # CRITICAL FIX: Restore the interpolation function to exactly match
                # the serial float math, preventing chaotic divergence!
                copyto!(F[n+1], local_sol(T[n+1]))
            end
        end

        # Save history if requested
        if saveiterates
            current_history_step = solution.iterates[k]
            for n = 1:N
                # CRITICAL FIX: Must deepcopy because local_sol mutates next iteration
                current_history_step[n] = deepcopy(solution[n])
            end
        end

        # Update weights of error function
        update!(weights, U, F)

        # Check convergence
        errors[k] = ψ(cache, k, weights)
        if errors[k] ≤ ϵ
            resize!(errors, k)
            if saveiterates
                resize!(solution.iterates, k)
            end
            break
        end

        # Correction step (serial)
        for n = k:N-1
            chunkproblem = copy(problem, U[n], T[n], T[n+1])
            chunkcoarsesolution = coarsesolver(chunkproblem)
            
            v = chunkcoarsesolution(T[n+1])
            
            # CRITICAL FIX: Rebind the array instead of mutating in-place.
            # This prevents U_previous in the error function from being overwritten!
            U[n+1] = v + F[n+1] - G[n+1]
            G[n+1] = v
        end
    end

    return solution
end
