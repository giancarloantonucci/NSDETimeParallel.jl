function getchunks(problem, solver::TimeParallelSolver)
    @↓ (t₀, tN) ← tspan = problem
    @↓ P = solver
    T = abs(tN - t₀)
    ΔT = T / P
    t₀:ΔT:tN
end

function getchunks(iterate::TimeParallelIterate)
    t₀ = iterate[1].t[1]
    tN = iterate[end].t[end]
    T = abs(tN - t₀)
    P = length(iterate)
    ΔT = T / P
    t₀:ΔT:tN
end

getchunks(solution::TimeParallelSolution) = getchunks(solution[end])
