"""
    coarseguess!(cache::PararealCache, solution::PararealSolution, problem::AbstractInitialValueProblem, parareal::Parareal)

computes a coarse initial guess and stores the result into `cache.U`.
"""
function coarseguess!(cache::PararealCache, solution::PararealSolution, problem::AbstractInitialValueProblem, parareal::Parareal)
    @↓ skips = cache
    @↓ U, T = solution.lastiterate
    @↓ u0, (t0, tN) ← tspan = problem
    @↓ coarsesolver = parareal
    @↓ N = parareal.parameters
    # stable sum
    for n = 1:N+1
        T[n] = (N - n + 1) / N * t0 + (n - 1) / N * tN
    end
    U[1] = u0
    for n = 1:N-1
        chunkproblem = copy(problem, U[n], T[n], T[n+1])
        chunkcoarsesolution = coarsesolver(chunkproblem)
        U[n+1] = chunkcoarsesolution(T[n+1])
        skips[n+1] = true
    end
    return cache
end

"""
    coarseguess(solution::PararealSolution, problem::AbstractInitialValueProblem, parareal::Parareal)

computes the coarse solution of `problem` and stores the result into a [`PararealCache`](@ref).
"""
function coarseguess(solution::PararealSolution, problem::InitialValueProblem, parareal::Parareal)
    cache = PararealCache(problem, parareal)
    coarseguess!(cache, solution, problem, parareal)
    return cache
end
