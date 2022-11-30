# function NSDEBase.cost(problem::AbstractInitialValueProblem, solver::AbstractTimeParallelSolver, solution::AbstractTimeParallelSolution)
#     @↓ (t0, tN) ← tspan = problem
#     @↓ finesolver, coarsolver, N = solver
#     @↓ hF ← h = finesolver.stepsize
#     @↓ hG ← h = coarsolver.stepsize
#     ξ = hG / hF
#     K = length(solution.errors)
#     l = steplength(solution.lastiterate)
#     return trunc(Int, l * (K + 1) / 2 * (N + ξ) / (N*ξ))
# end
