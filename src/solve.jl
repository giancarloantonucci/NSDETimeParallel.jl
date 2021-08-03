"""
    solve!(solution::TimeParallelSolution, problem, solver::TimeParallelSolver; kwargs...) :: TimeParallelSolution

returns the [`TimeParallelSolution`](@ref) of a problem, e.g. an [`InitialValueProblem`](@ref).
"""
NSDEBase.solve!(solution::TimeParallelSolution, problem, solver::TimeParallelSolver; kwargs...) = solver(solution, problem; kwargs...)

"""
    solve(problem, solver::TimeParallelSolver; kwargs...) :: TimeParallelSolution

returns the [`TimeParallelSolution`](@ref) of a problem, e.g. an [`InitialValueProblem`](@ref).
"""
NSDEBase.solve(problem, solver::TimeParallelSolver; kwargs...) = solver(problem; kwargs...)
