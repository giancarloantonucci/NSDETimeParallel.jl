"Serial implementation of Parareal (In-Memory)."
function parareal_serial!(
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

    # Initialization
    F[1] = G[1] = U[1]
    
    # Coarse run (serial)
    for n = 1:N-1
        chunkproblem = copy(problem, U[n], T[n], T[n+1])
        chunkcoarsesolution = coarsesolver(chunkproblem)
        G[n+1] = chunkcoarsesolution(T[n+1])
        if makeGs[n+1]
            U[n+1] = copy(G[n+1])
        end
    end

    # Main loop
    for k = 1:K
        # Fine run (Serial loop)
        for n = k:N
            chunkproblem = copy(problem, U[n], T[n], T[n+1])
            chunkfinesolution = finesolver(chunkproblem)
            
            # DIRECT STORAGE: Write directly into the solution object
            # This works because solution.lastiterate is a vector of chunks in memory
            solution[n] = chunkfinesolution 

            # Update F
            if n < N
                F[n+1] = chunkfinesolution(T[n+1])
            end
        end

        # Save history if requested
        if saveiterates
            # We copy the pointers from the current `lastiterate` (solution)
            # to the history vector `iterates[k]`
            current_history_step = solution.iterates[k]
            for n = 1:N
                current_history_step[n] = solution[n]
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
            U[n+1] = v + F[n+1] - G[n+1]
            G[n+1] = v
        end
    end

    return solution
end
