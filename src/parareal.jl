"""
    Parareal <: TimeParallelSolver

A composite type for the parareal algorithm.

# Constructors
```julia
Parareal(ℱ, 𝒢; P = 10, K = P, φ = ErrorCheck())
```

# Arguments
- `ℱ :: Union{Function, InitialValueSolver}` : fine solver.
- `𝒢 :: Union{Function, InitialValueSolver}` : coarse solver.
- `P :: Integer` : number of time chunks.
- `K :: Integer` : maximum number of iterations.
- `𝜑 :: Function` : error control function.
- `ϵ :: Real` : tolerance.

# Functions
- [`show`](@ref) : shows name and contents.
- [`summary`](@ref) : shows name.
"""
struct Parareal{error_check_T, ℱ_T, 𝒢_T, P_T, K_T} <: TimeParallelSolver
    error_check::error_check_T
    ℱ::ℱ_T
    𝒢::𝒢_T
    P::P_T
    K::K_T
end

function Parareal(ℱ::Function, 𝒢::Function; P = 10, K = P, φ = ErrorCheck())
    return Parareal(φ, ℱ, 𝒢, P, K)
end

function Parareal(finesolver::InitialValueSolver, coarsolver::InitialValueSolver; P = 10, K = P, φ = ErrorCheck())
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
    return Parareal(ℱ, 𝒢; P=P, K=K, φ=φ)
end

# ---------------------------------------------------------------------------- #
#                                   Functions                                  #
# ---------------------------------------------------------------------------- #

"""
    show(io::IO, parareal::Parareal)

prints a full description of `parareal` and its contents to a stream `io`.
"""
Base.show(io::IO, parareal::Parareal) = NSDEBase._show(io, parareal)
# function Base.show(io::IO, solver::Parareal)
#     print(io, "TimeParallelSolver:")
#     pad = get(io, :pad, "")
#     names = propertynames(solver)
#     N = length(names)
#     for (n, name) in enumerate(names)
#         field = getproperty(solver, name)
#         if string(name) == "ℱ" || string(name) == "𝒢"
#             print(io, "\n", pad, "   ‣ " * string(name) * " ≔ ")
#             summary(io, field)
#         else
#             if field !== nothing
#                 print(io, "\n", pad, "   ‣ " * string(name) * " ≔ ")
#                 show(IOContext(io, :pad => string(pad, "   ")), field)
#             end
#         end
#     end
# end

"""
    summary(io::IO, parareal::Parareal)

prints a brief description of `parareal` to a stream `io`.
"""
Base.summary(io::IO, parareal::Parareal) = NSDEBase._summary(io, parareal)

# ---------------------------------------------------------------------------- #
#                                    Methods                                   #
# ---------------------------------------------------------------------------- #

"""
    coarseguess!(solution::TimeParallelSolution, problem, u0, t0, tN, solver::Parareal)

computes the coarse solution of a `problem`, e.g. an [`InitialValueProblem`](@ref), as part of the first iteration of [`Parareal`](@ref).
"""
function coarseguess!(solution::TimeParallelSolution, problem, u0, t0, tN, solver::Parareal)
    @↓ 𝒢, P = solver
    @↓ U, T = solution
    T[1] = t0
    for n in 1:P
        # more stable sum
        T[n+1] = (1 - n / P) * t0 + n * tN / P
    end
    U[1] = u0
    for n = 1:P
        # subproblem = IVP(rhs, U[n], T[n], T[n+1])
        # chunk = solve(problem, coarsolver)
        chunk = 𝒢(problem, U[n], T[n], T[n+1])
        U[n+1] = chunk.u[end]
    end
    @↑ solution = U, T
end

"""
    coarseguess!(solution::TimeParallelSolution, problem, solver::Parareal)

computes the coarse solution of a `problem`, e.g. an [`InitialValueProblem`](@ref), as part of the first iteration of [`Parareal`](@ref).
"""
function coarseguess!(solution::TimeParallelSolution, problem, solver::Parareal)
    @↓ u0, (t0, tN) ← tspan = problem
    coarseguess!(solution, problem, u0, t0, tN, solver)
end

function parareal_serial!(solution::TimeParallelSolution, problem, solver::Parareal)
    @↓ iterates, φ, U, T = solution
    @↓ ℱ, 𝒢, P, K = solver
    @↓ 𝜑, ϵ, Λ, updateΛ = solver.error_check
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
        # @↑ solution[k] = U .← U
        solution[k].U .= U
        for n = 1:k-1
            solution[k][n] = solution[k-1][n]
        end
        # fine run (parallelisable)
        for n = k:P
            chunk = ℱ(problem, U[n], T[n], T[n+1])
            solution[k][n] = chunk
            F[n+1] = chunk.u[end]
        end
        solution[k].F .= F
        # update Lipschitz constant
        Λ = updateΛ ? update_Lipschitz(Λ, U, F) : Λ
        # check convergence
        φ[k] = 𝜑(solution, k, Λ)
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

function parareal_distributed!(solution::TimeParallelSolution, problem, solver::Parareal)
    @↓ iterates, φ, U, T = solution
    @↓ ℱ, 𝒢, P, K = solver
    @↓ 𝜑, ϵ, Λ, updateΛ = solver.error_check
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
        # @↑ solution[k] = U .← U
        solution[k].U .= U
        for n = 1:k-1
            solution[k][n] = solution[k-1][n]
        end
        # fine run (with Julia's Distributed.jl)
        @sync for n = k:P
            @async F[n+1] = remotecall_fetch(getF, n, problem, U[n], T[n], T[n+1])
        end
        solution[k].F .= F
        # update Lipschitz constant
        Λ = updateΛ ? update_Lipschitz(Λ, U, F) : Λ
        # check convergence
        φ[k] = 𝜑(solution, k, Λ)
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

function (solver::Parareal)(solution::TimeParallelSolution, problem; mode = "SERIAL")
    if nprocs() == 1 || mode == "SERIAL"
        parareal_serial!(solution, problem, solver)
    elseif nprocs() > 1 && mode == "DISTRIBUTED"
        parareal_distributed!(solution, problem, solver)
    elseif nprocs() > 1 && mode == "MPI"
        # parareal_mpi!(solution, problem, solver)
    end
    return solution
end

function (solver::Parareal)(problem; mode = "SERIAL")
    solution = TimeParallelSolution(problem, solver)
    coarseguess!(solution, problem, solver)
    solver(solution, problem; mode=mode)
    return solution
end
