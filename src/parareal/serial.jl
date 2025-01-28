"Serial implementation of Parareal."
function parareal_serial!(cache::PararealCache, solution::PararealSolution, problem::AbstractInitialValueProblem, parareal::Parareal)
    @↓ skips, F, G = cache
    @↓ errors, iterates = solution
    @↓ U, T = solution.lastiterate
    @↓ u0, (t0, tN) ← tspan = problem
    @↓ finesolver, coarsesolver, saveiterates = parareal
    @↓ N, K = parareal.parameters
    @↓ weights, ψ, ϵ = parareal.tolerance

    # Initialization
    for n = 1:N+1
        T[n] = (N - n + 1) / N * t0 + (n - 1) / N * tN # stable sum
    end
    F[1] = G[1] = U[1] = u0

    # Coarse run (serial)
    for n = 1:N-1
        chunkproblem = copy(problem, U[n], T[n], T[n+1])
        chunkcoarsesolution = coarsesolver(chunkproblem)
        G[n+1] = chunkcoarsesolution(T[n+1])
        if !skips[n+1]
            U[n+1] = G[n+1]
        end
    end

    # Main loop
    for k in 1:K

        # Fine run (parallelisable)
        for n = k:N
            chunkproblem = copy(problem, U[n], T[n], T[n+1])
            chunkfinesolution = finesolver(chunkproblem)
            solution[n] = chunkfinesolution
            if n < N
                F[n+1] = chunkfinesolution(T[n+1])
            end
        end

        # Save iterates
        if saveiterates
            iterates[k].U .= U
            iterates[k].T .= T
            for n = 1:k-1
                iterates[k][n] = iterates[k-1][n]
            end
            for n = k:N
                iterates[k][n] = solution[n]
            end
        end

        # Update weights of error function
        update!(weights, U, F)

        # Check convergence
        errors[k] = ψ(cache, solution, k, weights)
        if errors[k] ≤ ϵ
            resize!(errors, k) # self-updates inside solution
            if saveiterates
                resize!(iterates, k) # self-updates inside solution
            end
            break
        end

        # correction step (serial)
        for n = k:N-1
            chunkproblem = copy(problem, U[n], T[n], T[n+1])
            chunkfinesolution = coarsesolver(chunkproblem)
            v = chunkfinesolution(T[n+1])
            U[n+1] = v + F[n+1] - G[n+1]
            G[n+1] = v
        end
        
    end

    return solution
end
