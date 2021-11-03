function parareal_serial!(cache::PararealCache, solution::PararealSolution, problem::AbstractInitialValueProblem, parareal::Parareal)
    @↓ U, T, F, G = cache
    @↓ errors, saveiterates = solution
    @↓ u0, (t0, tN) ← tspan = problem
    @↓ finesolver, coarsolver, P, K = parareal
    @↓ weights, ψ, ϵ = parareal.control
    T[1] = t0
    for n = 1:P
        # stable sum
        T[n+1] = (1-n) * t0 / P + n * tN / P
    end
    F[1] = G[1] = U[1] = u0
    # coarse guess
    for n = 1:P
        chunkproblem = subproblemof(problem, U[n], T[n], T[n+1])
        chunksolution = coarsolver(chunkproblem)
        U[n+1] = chunksolution.u[end]
    end
    G .= U
    for k in 1:K
        if saveiterates
            for n = 1:(k - 1)
                solution[k][n] = solution[k-1][n]
            end
        end
        # fine run (parallelisable)
        for n = k:P
            chunkproblem = subproblemof(problem, U[n], T[n], T[n+1])
            chunksolution = finesolver(chunkproblem)
            F[n + 1] = chunksolution.u[end]
            if saveiterates
                solution[k][n] = chunksolution
            else
                solution[n] = chunksolution
            end
        end
        # check convergence
        update!(weights, U, F)
        errors[k] = ψ(cache, solution, k, weights)
        if errors[k] ≤ ϵ
            resize!(errors, k)
            resize!(solution, k)
            break
        end
        # serial update
        for n = k:P
            chunkproblem = subproblemof(problem, U[n], T[n], T[n+1])
            chunksolution = coarsolver(chunkproblem)
            U[n+1] = chunksolution.u[end] + F[n+1] - G[n+1]
            G[n+1] = chunksolution.u[end]
        end
    end
    return solution
end
