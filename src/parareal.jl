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

include("parareal/coarse.jl")
include("parareal/serial.jl")
include("parareal/distributed.jl")
include("parareal/mpi.jl")

function (solver::Parareal)(solution::TimeParallelSolution, problem; mode = "SERIAL")
    if nprocs() == 1 || mode == "SERIAL"
        parareal_serial!(solution, problem, solver)
    elseif nprocs() > 1 && mode == "DISTRIBUTED"
        parareal_distributed!(solution, problem, solver)
    elseif nprocs() > 1 && mode == "MPI"
        parareal_mpi!(solution, problem, solver)
    end
    return solution
end

function (solver::Parareal)(problem; mode = "SERIAL")
    solution = TimeParallelSolution(problem, solver)
    coarseguess!(solution, problem, solver)
    solver(solution, problem; mode=mode)
    return solution
end
