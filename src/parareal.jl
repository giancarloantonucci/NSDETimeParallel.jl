struct Parareal{objective_T, mode_T, ℱ_T, 𝒢_T, P_T, K_T} <: TimeParallelSolver
    objective::objective_T
    mode::mode_T
    ℱ::ℱ_T
    𝒢::𝒢_T
    P::P_T
    K::K_T
end

function Parareal(ℱ::Function, 𝒢::Function; 𝜑 = 𝜑₁, ϵ = 1e-12, P = 10, K = P, mode = "SERIAL")
    objective = ErrorFunction(𝜑, ϵ)
    return Parareal(objective, mode, ℱ, 𝒢, P, K)
end

function Parareal(finesolver::InitialValueSolver, coarsolver::InitialValueSolver; 𝜑 = 𝜑₁, ϵ = 1e-12, P = 10, K = P, mode = "SERIAL")
    objective = ErrorFunction(𝜑, ϵ)
    function ℱ(problem, u0, tspan)
        subproblem = IVP(problem.rhs, u0, tspan)
        solve(subproblem, finesolver)
    end
    ℱ(problem, u0, t0, tN) = ℱ(problem, u0, (t0, tN))
    function 𝒢(problem, u0, tspan)
        subproblem = IVP(problem.rhs, u0, tspan)
        solve(subproblem, coarsolver)
    end
    𝒢(problem, u0, t0, tN) = 𝒢(problem, u0, (t0, tN))
    @everywhere begin
        finesolver = $finesolver
        coarsolver = $coarsolver
        function ℱ(problem, u0, tspan)
            subproblem = IVP(problem.rhs, u0, tspan)
            solve(subproblem, finesolver)
        end
        ℱ(problem, u0, t0, tN) = ℱ(problem, u0, (t0, tN))
        function 𝒢(problem, u0, tspan)
            subproblem = IVP(problem.rhs, u0, tspan)
            solve(subproblem, coarsolver)
        end
        𝒢(problem, u0, t0, tN) = 𝒢(problem, u0, (t0, tN))
    end
    return Parareal(objective, mode, ℱ, 𝒢, P, K)
end

function solve!(solution::TimeParallelSolution, problem, solver::Parareal)
    @↓ iterates, φ, U, T = solution
    @↓ mode, ℱ, 𝒢, P, K = solver
    @↓ 𝜑, ϵ = solver.objective
    # coarse guess
    G = similar(U); G[1] = U[1]
    for n = 1:P
        chunk = 𝒢(problem, U[n], T[n], T[n+1])
        G[n+1] = chunk.u[end]
    end
    # main loop
    F = similar(U); F[1] = U[1]
    getF(args...) = ℱ(args...).u[end]
    for k = 1:K
        # fine run (parallelisable)
        if mode == "SERIAL" || nprocs() == 1
            for n = k:P
                chunk = ℱ(problem, U[n], T[n], T[n+1])
                solution[k][n] = chunk
                F[n+1] = chunk.u[end]
            end
        elseif mode == "PARALLEL" && nprocs() > 1
            @sync for n = k:P
                @async F[n+1] = remotecall_fetch(getF, n, problem, U[n], T[n], T[n+1])
            end
        end
        # check convergence
        φ[k] = 𝜑(U, F, T)
        if φ[k] ≤ ϵ
            resize!(iterates, k)
            resize!(φ, k)
            break
        end
        # update (serial)
        for n = k:P
            chunk = 𝒢(problem, U[n], T[n], T[n+1])
            U[n+1] = chunk.u[end] + F[n+1] - G[n+1]
            G[n+1] = chunk.u[end]
        end
        @↑ solution = U, T
    end
    return solution
end
