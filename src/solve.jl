"""
    coarseguess!(solution::TimeParallelSolution, problem, solver::TimeParallelSolver)
    coarseguess!(solution::TimeParallelSolution, problem, u0, t0, tN, solver::TimeParallelSolver)

computes the coarse solution of a `problem`, e.g. an [`InitialValueProblem`](@ref), for the first iteration of a [`TimeParallelSolver`](@ref)
"""
function coarseguess!(solution::TimeParallelSolution, problem, u0, t0, tN, solver::TimeParallelSolver)
    @↓ 𝒢, P = solver
    @↓ U, T = solution
    T[1] = t0
    for n in 1:P
        # more stable sum
        T[n+1] = (1 - n / P) * t0 + n * tN / P
    end
    U[1] = u0
    for n = 1:P
        chunk = 𝒢(problem, U[n], T[n], T[n+1])
        U[n+1] = chunk.u[end]
    end
    @↑ solution = U, T
end

function coarseguess!(solution::TimeParallelSolution, problem, solver::TimeParallelSolver)
    @↓ u0, (t0, tN) ← tspan = problem
    coarseguess!(solution, problem, u0, t0, tN, solver)
end

"""
    solve!(solution::TimeParallelSolution, problem, solver::TimeParallelSolver; mode::String = "SERIAL")

returns the `TimeParallelSolution` of a problem, e.g. an [`InitialValueProblem`](@ref).
"""
function NSDEBase.solve!(solution::TimeParallelSolution, problem, solver::TimeParallelSolver; mode = "SERIAL")
    coarseguess!(solution, problem, solver)
    if nprocs() == 1 || mode == "SERIAL"
        solve_serial!(solution, problem, solver)
    elseif nprocs() > 1 && mode == "DISTRIBUTED"
        solve_distributed!(solution, problem, solver)
    elseif nprocs() > 1 && mode == "MPI"
        # solve_mpi!(solution, problem, solver)
    end
    return solution
end

"""
    solve(problem, solver::TimeParallelSolver; mode::String = "SERIAL")

returns the `TimeParallelSolution` of a problem, e.g. an [`InitialValueProblem`](@ref).
"""
function NSDEBase.solve(problem, solver::TimeParallelSolver; mode = "SERIAL")
    solution = TimeParallelSolution(problem, solver)
    solve!(solution, problem, solver; mode=mode)
    return solution
end
