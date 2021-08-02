"""
    Parareal{error_check_T, ℱ_T, 𝒢_T, P_T, K_T} <: TimeParallelSolver

returns a constructor for the [`TimeParallelSolver`](@ref) based on the parareal algorithm.

---

    Parareal(ℱ::Function, 𝒢::Function; P = 10, K = P, 𝜑 = 𝜑₁, ϵ = 1e-12)

returns a [`Parareal`](@ref) with:
- `ℱ :: Function` : fine solver.
- `𝒢 :: Function` : coarse solver.
- `P :: Integer`  : number of time chunks.
- `K :: Integer`  : maximum number of iterations.
- `𝜑 :: Function` : error control function.
- `ϵ :: Real`     : tolerance.

---

    Parareal(finesolver::InitialValueSolver, coarsolver::InitialValueSolver; P = 10, K = P, 𝜑 = 𝜑₁, ϵ = 1e-12)

returns a [`Parareal`](@ref) from a `finesolver` and a `coarsolver`.
"""
struct Parareal{error_check_T, ℱ_T, 𝒢_T, P_T, K_T} <: TimeParallelSolver
    error_check::error_check_T
    ℱ::ℱ_T
    𝒢::𝒢_T
    P::P_T
    K::K_T
end

function Parareal(ℱ::Function, 𝒢::Function; P = 10, K = P, 𝜑 = 𝜑₁, ϵ = 1e-12)
    error_check = ErrorCheck(𝜑, ϵ)
    return Parareal(error_check, ℱ, 𝒢, P, K)
end

function Parareal(finesolver::InitialValueSolver, coarsolver::InitialValueSolver; P = 10, K = P, 𝜑 = 𝜑₁, ϵ = 1e-12)
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
    return Parareal(ℱ, 𝒢; P=P, K=K, 𝜑=𝜑, ϵ=ϵ)
end

function solve_serial!(solution::TimeParallelSolution, problem, solver::Parareal)
    @↓ iterates, φ, U, T = solution
    @↓ ℱ, 𝒢, P, K = solver
    @↓ 𝜑, ϵ = solver.error_check
    # coarse guess
    G = similar(U)
    G[1] = U[1]
    for n = 1:P
        chunk = 𝒢(problem, U[n], T[n], T[n+1])
        G[n+1] = chunk.u[end]
    end
    # main loop
    F = similar(U)
    F[1] = U[1]
    for k = 1:K
        for n = 1:k-1
            solution[k][n] = solution[k-1][n]
        end
        # fine run (parallelisable)
        for n = k:P
            chunk = ℱ(problem, U[n], T[n], T[n+1])
            solution[k][n] = chunk
            F[n+1] = chunk.u[end]
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

function solve_distributed!(solution::TimeParallelSolution, problem, solver::Parareal)
    @↓ iterates, φ, U, T = solution
    @↓ ℱ, 𝒢, P, K = solver
    @↓ 𝜑, ϵ = solver.error_check
    # coarse guess
    G = similar(U)
    G[1] = U[1]
    for n = 1:P
        chunk = 𝒢(problem, U[n], T[n], T[n+1])
        G[n+1] = chunk.u[end]
    end
    # main loop
    F = similar(U)
    F[1] = U[1]
    getF(args...) = ℱ(args...).u[end]
    for k = 1:K
        # fine run (parallelisable)
        @sync for n = k:P
            @async F[n+1] = remotecall_fetch(getF, n, problem, U[n], T[n], T[n+1])
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
