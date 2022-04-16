function fulllength(iterate::AbstractTimeParallelIterate)
    l = 0
    for chunk in iterate.chunks
        l += length(chunk)
    end
    return l-1
end

# function NSDEBase.cost(problem::AbstractInitialValueProblem, solver::AbstractTimeParallelSolver, solution::AbstractTimeParallelSolution)
#     @↓ (t0, tN) ← tspan = problem
#     @↓ finesolver, coarsolver, P = solver
#     @↓ hF ← h = finesolver.stepsize
#     @↓ hG ← h = coarsolver.stepsize
#     ξ = hG / hF
#     K = length(solution.errors)
#     l = steplength(solution.lastiterate)
#     return trunc(Int, l * (K + 1) / 2 * (P + ξ) / (P*ξ))
# end
