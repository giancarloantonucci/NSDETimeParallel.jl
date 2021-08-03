"""
    getchunks(problem, solver::TimeParallelSolver) :: UnitRange

returns the time-chunk nodes.
"""
function getchunks(problem, solver::TimeParallelSolver)
    @↓ (t0, tN) ← tspan = problem
    @↓ P = solver
    T = abs(tN - t0)
    ΔT = T / P
    t0:ΔT:tN
end

"""
    getchunks(iterate::TimeParallelIterate) :: UnitRange

returns the time-chunk nodes.
"""
function getchunks(iterate::TimeParallelIterate)
    t0 = iterate[1].t[1]
    tN = iterate[end].t[end]
    T = abs(tN - t0)
    P = length(iterate)
    ΔT = T / P
    t0:ΔT:tN
end

"""
    getchunks(solution::TimeParallelSolution) :: UnitRange

returns the time-chunk nodes.
"""
getchunks(solution::TimeParallelSolution) = getchunks(solution[end])
