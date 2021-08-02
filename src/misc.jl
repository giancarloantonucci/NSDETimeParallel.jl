"""
    getchunks(problem, solver::TimeParallelSolver)
    getchunks(iterate::TimeParallelIterate)
    getchunks(solution::TimeParallelSolution)

returns a `UnitRange` of the nodes at which the time domain has been divided into chunks.
"""
function getchunks(problem, solver::TimeParallelSolver)
    @↓ (t0, tN) ← tspan = problem
    @↓ P = solver
    T = abs(tN - t0)
    ΔT = T / P
    t0:ΔT:tN
end

function getchunks(iterate::TimeParallelIterate)
    t0 = iterate[1].t[1]
    tN = iterate[end].t[end]
    T = abs(tN - t0)
    P = length(iterate)
    ΔT = T / P
    t0:ΔT:tN
end

getchunks(solution::TimeParallelSolution) = getchunks(solution[end])
