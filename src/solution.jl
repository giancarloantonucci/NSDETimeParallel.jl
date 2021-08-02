"""
    TimeParallelSolution{iterates_T, φ_T, U_T, T_T} <: InitialValueSolution

returns a constructor for the numerical solution of an [`InitialValueProblem`](@ref) obtained with a [`TimeParallelSolver`](@ref).

---

    TimeParallelSolution(problem, solver::TimeParallelSolver)

returns an initialised [`TimeParallelSolution`](@ref) given a `problem`, e.g. an [`InitialValueProblem`](@ref), and a [`TimeParallelSolver`](@ref).
"""
mutable struct TimeParallelSolution{iterates_T, φ_T, U_T, T_T} <: InitialValueSolution
    iterates::iterates_T
    φ::φ_T
    U::U_T
    T::T_T
end

function TimeParallelSolution(problem, solver::TimeParallelSolver)
    @↓ u0, tspan = problem
    @↓ P = solver
    @↓ ϵ = solver.error_check
    iterates = [TimeParallelIterate(problem, solver) for i = 1:P]
    φ = Vector{typeof(ϵ)}(undef, P)
    U = Vector{typeof(u0)}(undef, P+1)
    T = Vector{eltype(tspan)}(undef, P+1)
    return TimeParallelSolution(iterates, φ, U, T)
end

# solution[k] ≡ solution.iterates[k]
Base.length(solution::TimeParallelSolution) = length(solution.iterates)
Base.getindex(solution::TimeParallelSolution, k::Int) = solution.iterates[k]
Base.setindex!(solution::TimeParallelSolution, value::TimeParallelIterate, k::Int) = solution.iterates[k] = value
Base.lastindex(solution::TimeParallelSolution) = lastindex(solution.iterates)

(solution::TimeParallelSolution)(t::Real) = solution[end](t)
(solution::TimeParallelSolution)(t::Real, n::Integer) = solution[n](t)
