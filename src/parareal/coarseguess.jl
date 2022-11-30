"""
    coarseguess!(cache::PararealCache, problem::AbstractInitialValueProblem, parareal::Parareal)

computes a coarse initial guess and stores the result into `cache.U`.
"""
function coarseguess!(cache::PararealCache, problem::AbstractInitialValueProblem, parareal::Parareal)
    @↓ skips, U, T = cache
    @↓ u0, (t0, tN) ← tspan = problem
    @↓ coarsolver = parareal
    @↓ N = parareal.parameters
    for n = 1:N+1
        T[n] = (N - n + 1) / N * t0 + (n - 1) / N * tN # stable sum
    end
    U[1] = u0
    for n = 1:N-1
        chunkproblem = subproblemof(problem, U[n], T[n], T[n+1])
        chunksolution = coarsolver(chunkproblem)
        U[n+1] = chunksolution(T[n+1])
        skips[n+1] = true
    end
    return cache
end

"""
    coarseguess(problem::AbstractInitialValueProblem, parareal::Parareal)

computes the coarse solution of `problem` and stores the result into a [`PararealCache`](@ref).
"""
function coarseguess(problem::InitialValueProblem, parareal::Parareal)
    cache = PararealCache(problem, parareal)
    coarseguess!(cache, problem, parareal)
    return cache
end
