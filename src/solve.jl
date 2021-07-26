function coarseguess!(solution::TimeParallelSolution, problem, u0, t0, tN, solver::TimeParallelSolver)
    @↓ 𝒢, P = solver
    @↓ U, T = solution
    ΔT = (tN - t0) / P
    T[1] = t0
    for n = 1:P
        T[n+1] = T[n] + ΔT
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

function NSDEBase.solve(problem, solver::TimeParallelSolver; mode = "SERIAL")
    solution = TimeParallelSolution(problem, solver)
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
