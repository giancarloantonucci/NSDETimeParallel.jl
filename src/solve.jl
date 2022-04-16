"""
    solve!(solution::AbstractTimeParallelSolution, problem, solver::AbstractTimeParallelSolver; kwargs...) :: AbstractTimeParallelSolution

computes the `solution` of `problem` using `solver`.
"""
function NSDEBase.solve!(solution::AbstractTimeParallelSolution, problem::AbstractInitialValueProblem, solver::AbstractTimeParallelSolver; kwargs...)
    return solver(solution, problem; kwargs...)
end

"""
    solve(problem, solver::AbstractTimeParallelSolver; kwargs...) :: AbstractTimeParallelSolution

computes the solution of `problem` using `solver`.
"""
function NSDEBase.solve(problem::AbstractInitialValueProblem, solver::AbstractTimeParallelSolver; kwargs...)
    return solver(problem; kwargs...)
end
