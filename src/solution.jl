mutable struct TimeParallelIterate{chunks_T, U_T, T_T}
    chunks::chunks_T
    U::U_T
    T::T_T
end

function TimeParallelIterate(problem, solver::TimeParallelSolver)
    @↓ u0, t0 ← tspan[1] = problem
    @↓ P = solver
    chunks = Vector{Any}(undef, P)
    U = Vector{typeof(u0)}(undef, P+1)
    T = Vector{typeof(t0)}(undef, P+1)
    TimeParallelIterate(chunks, U, T)
end

# solution[k][n] ≡ solution.iterates[k].chunks[n]
Base.length(iterate::TimeParallelIterate) = length(iterate.chunks)
Base.getindex(iterate::TimeParallelIterate, n::Int) = iterate.chunks[n]
Base.setindex!(iterate::TimeParallelIterate, value, n::Int) = iterate.chunks[n] = value
Base.lastindex(iterate::TimeParallelIterate) = lastindex(iterate.chunks)

mutable struct TimeParallelSolution{iterates_T, φ_T, U_T, T_T} <: InitialValueSolution
    iterates::iterates_T
    φ::φ_T
    U::U_T
    T::T_T
end

function TimeParallelSolution(problem, solver::TimeParallelSolver)
    @↓ u0, tspan = problem
    @↓ P = solver
    @↓ ϵ = solver.objective
    iterates = fill(TimeParallelIterate(problem, solver), P)
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
