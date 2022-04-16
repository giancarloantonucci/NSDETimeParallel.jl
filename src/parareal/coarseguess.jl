"""
    coarseguess!(cache::PararealCache, problem::AbstractInitialValueProblem, parareal::Parareal)

computes the coarse solution of `problem` and stores the result into `cache`.
"""
function coarseguess!(cache::PararealCache, problem::AbstractInitialValueProblem, parareal::Parareal)
    @↓ U, G, T = cache
    @↓ u0, (t0, tN) ← tspan = problem
    @↓ coarsolver, P = parareal
    # stable sum
    for n = 1:P+1
        T[n] = (P - n + 1) / P * t0 + (n - 1) / P * tN
    end
    G[1] = U[1] = u0
    for n = 1:P-1
        chunkproblem = subproblemof(problem, U[n], T[n], T[n+1])
        G[n+1] = U[n+1] = coarsolver(chunkproblem).u[end]
    end
    return cache
end

function coarseguess(problem::InitialValueProblem, parareal::Parareal)
    cache = PararealCache(problem, parareal)
    coarseguess!(cache, problem, parareal)
    return cache
end
