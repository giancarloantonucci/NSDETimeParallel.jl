"""
    coarseguess!(solution::TimeParallelSolution, problem, u0, t0, tN, solver::Parareal)

computes the coarse solution of a `problem`, e.g. an [`InitialValueProblem`](@ref), as part of the first iteration of [`Parareal`](@ref).
"""
function coarseguess!(solution::TimeParallelSolution, problem, u0, t0, tN, solver::Parareal)
    @↓ 𝒢, P = solver
    @↓ U, T = solution
    T[1] = t0
    for n in 1:P
        # more stable sum
        T[n+1] = (1 - n / P) * t0 + n * tN / P
    end
    U[1] = u0
    for n = 1:P
        # subproblem = IVP(rhs, U[n], T[n], T[n+1])
        # chunk = solve(problem, coarsolver)
        chunk = 𝒢(problem, U[n], T[n], T[n+1])
        U[n+1] = chunk.u[end]
    end
    @↑ solution = U, T
end

"""
    coarseguess!(solution::TimeParallelSolution, problem, solver::Parareal)

computes the coarse solution of a `problem`, e.g. an [`InitialValueProblem`](@ref), as part of the first iteration of [`Parareal`](@ref).
"""
function coarseguess!(solution::TimeParallelSolution, problem, solver::Parareal)
    @↓ u0, (t0, tN) ← tspan = problem
    coarseguess!(solution, problem, u0, t0, tN, solver)
end
