mutable struct MovingWindowSolver{𝒫_T, τ_T, Δτ_T} <: InitialValueSolver
    𝒫::𝒫_T
    τ::τ_T
    Δτ::Δτ_T
end

# function MovingWindowSolver(parallelsolver::TimeParallelSolver, τ, Δτ)
#     if Δτ < τ / parallelsolver.P
#         error("Select Δτ ≥ τ / P!")
#     elseif Δτ > τ
#         error("Select Δτ ≤ τ!")
#     end
#     new(parallelsolver, τ, Δτ)
# end

MovingWindowSolver(𝒫; τ, Δτ) = MovingWindowSolver(𝒫, τ, Δτ)
@doc (@doc MovingWindowSolver) MoWA(args...; kwargs...) = MovingWindowSolver(args...; kwargs...)

mutable struct MovingWindowSolution{windows_T} <: InitialValueSolution
    windows::windows_T
end

function MovingWindowSolution(problem, solver::MovingWindowSolver)
    @↓ (t0, tN) ← tspan = problem
    @↓ τ, Δτ = solver
    M = 1 # if solver.Δτ == 0 then T = τ and M = 1
    if Δτ > 0
        M = trunc(Int, (tN - τ) / Δτ) + 1
    end
    windows = Vector{TimeParallelSolution}(undef, M)
    return MovingWindowSolution(windows)
end

Base.length(solution::MovingWindowSolution) = length(solution.windows)
Base.getindex(solution::MovingWindowSolution, m::Int) = solution.windows[m] # solution[m] ≡ solution.windows[m]
Base.setindex!(solution::MovingWindowSolution, value::TimeParallelSolution, m::Int) = solution.windows[m] = value
Base.lastindex(solution::MovingWindowSolution) = lastindex(solution.windows)

function solve!(solution::MovingWindowSolution, problem, solver::MovingWindowSolver)
    @↓ u0, (t0, tN) ← tspan = problem
    @↓ 𝒫, τ, Δτ = solver
    @↓ 𝒢, P = 𝒫
    for m = 1:length(solution)
        solution[m] = TimeParallelSolution(problem, 𝒫)
        tmp = solution[m]
        @↓ U, T = tmp
        if m == 1
            coarseguess!(solution[m], problem, u0, t0, t0 + τ, 𝒫)
        else
            ΔP = trunc(Int, P * Δτ / τ)
            N = P - ΔP + 1
            for n = 1:length(T)
                T[n] = solution[m-1].T[n] + Δτ
            end
            for n = 1:N
                U[n] = solution[m-1].U[ΔP+n]
            end
            for n = N:P
                chunk = 𝒢(problem, U[n], T[n], T[n+1])
                U[n+1] = chunk.u[end]
            end
        end
        solve_serial!(solution[m], problem, 𝒫)
    end
    solution
end

function NSDEBase.solve(problem, solver::MovingWindowSolver)
    solution = MovingWindowSolution(problem, solver)
    solve!(solution, problem, solver)
    solution
end
