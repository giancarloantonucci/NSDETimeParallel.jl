"""
    coarseguess!(cache::PararealCache, problem::AbstractInitialValueProblem, parareal::Parareal)

computes the coarse solution of an [`AbstractInitialValueProblem`](@ref), and saves the outcome into a [`PararealCache`](@ref).
"""
function coarseguess!(cache::PararealCache, problem::AbstractInitialValueProblem, parareal::Parareal)
    @↓ U, T, G = cache
    @↓ u0, (t0, tN) ← tspan = problem
    @↓ coarsolver, P = parareal
    T[1] = t0
    for n = 1:P
        # stable sum
        T[n + 1] = (1 - n) * t0 / P + n * tN / P
    end
    G[1] = u0
    for n = 1:P
        chunkproblem = subproblemof(problem, G[n], T[n], T[n+1])
        chunksolution = coarsolver(chunkproblem)
        G[n + 1] = chunksolution.u[end]
    end
    U .= G
    return cache
end
