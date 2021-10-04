"""
    TimeParallelSolution <: InitialValueSolution

A composite type for the solution of an [`InitialValueProblem`](@ref) obtained with a [`TimeParallelSolver`](@ref).

# Constructors
```julia
TimeParallelSolution(problem, solver)
```

# Arguments
- `problem` : initial value problem, e.g. an [`InitialValueProblem`](@ref).
- `solver :: TimeParallelSolver`.

# Functions
- [`getindex`](@ref) : get iterate.
- [`lastindex`](@ref) : last index.
- [`length`](@ref) : number of iterates.
- [`setindex!`](@ref) : set iterate.
- [`show`](@ref) : shows name and contents.
- [`summary`](@ref) : shows name.
"""
mutable struct TimeParallelSolution{iterates_T, φ_T, U_T, T_T} <: InitialValueSolution
    iterates::iterates_T
    φ::φ_T
    U::U_T
    T::T_T
end

function TimeParallelSolution(problem, solver::TimeParallelSolver)
    @↓ u0, tspan = problem
    @↓ P, K = solver
    @↓ ϵ = solver.error_check
    iterates = [TimeParallelIterate(problem, solver) for i = 1:K]
    φ = Vector{typeof(ϵ)}(undef, P)
    U = Vector{typeof(u0)}(undef, P+1)
    T = Vector{eltype(tspan)}(undef, P+1)
    return TimeParallelSolution(iterates, φ, U, T)
end

# ---------------------------------------------------------------------------- #
#                                   Functions                                  #
# ---------------------------------------------------------------------------- #

# solution[k] ≡ solution.iterates[k]

"""
    length(solution::TimeParallelSolution)

returns the number of iterates of `solution`.
"""
Base.length(solution::TimeParallelSolution) = length(solution.iterates)

"""
    getindex(solution::TimeParallelSolution, k::Integer)

returns the `n`-th iterate of a [`TimeParallelSolution`](@ref).
"""
Base.getindex(solution::TimeParallelSolution, k::Integer) = solution.iterates[k]

"""
    setindex!(solution::TimeParallelSolution, iterate::TimeParallelIterate, k::Integer)

stores a [`TimeParallelIterate`](@ref) into the `n`-th iterate of a [`TimeParallelSolution`](@ref).
"""
Base.setindex!(solution::TimeParallelSolution, iterate::TimeParallelIterate, k::Integer) = solution.iterates[k] = iterate

"""
    lastindex(solution::TimeParallelSolution)

returns the last index of `solution`.
"""
Base.lastindex(solution::TimeParallelSolution) = lastindex(solution.iterates)

"""
    show(io::IO, solution::TimeParallelSolution)

prints a full description of `solution` and its contents to a stream `io`.
"""
Base.show(io::IO, solution::TimeParallelSolution) = NSDEBase._show(io, solution)

"""
    summary(io::IO, solution::TimeParallelSolution)

prints a brief description of `solution` to a stream `io`.
"""
Base.summary(io::IO, solution::TimeParallelSolution) = NSDEBase._summary(io, solution)

function RungeKutta.extract(iterate::TimeParallelIterate, i::Integer)
    N = length(iterate)
    @↓ u1 ← u, t1 ← t = RungeKutta.extract(iterate[1], i)
    u = eltype(u1)[]
    t = eltype(t1)[]
    append!(u, u1[1:end])
    append!(t, t1[1:end])
    for n = 2:N
        @↓ u1 ← u, t1 ← t = RungeKutta.extract(iterate[n], i)
        append!(u, u1[2:end])
        append!(t, t1[2:end])
    end
    return RungeKutta.RungeKuttaSolution(u, t)
end

RungeKutta.extract(solution::TimeParallelSolution, i::Integer) = RungeKutta.extract(solution[end], i)

# ---------------------------------------------------------------------------- #
#                                    Methods                                   #
# ---------------------------------------------------------------------------- #

(solution::TimeParallelSolution)(t::Real) = solution[end](t)
(solution::TimeParallelSolution)(t::Real, n::Integer) = solution[n](t)
