"""
    TimeParallelIterate

A composite type for the single iterations within a [`TimeParallelSolution`](@ref).

# Constructors
```julia
TimeParallelIterate(problem, solver)
```

# Arguments
- `problem` : initial value problem, e.g. an [`InitialValueProblem`](@ref).
- `solver :: TimeParallelSolver`.

# Functions
- [`getindex`](@ref) : get chunk.
- [`lastindex`](@ref) : last index.
- [`length`](@ref) : number of chunks.
- [`setindex!`](@ref) : set chunk.
- [`show`](@ref) : shows name and contents.
- [`summary`](@ref) : shows name.
"""
mutable struct TimeParallelIterate{chunks_T, U_T, T_T, F_T}
    chunks::chunks_T
    U::U_T
    T::T_T
    F::F_T
end

function TimeParallelIterate(problem, solver::TimeParallelSolver)
    @↓ u0, tspan = problem
    @↓ P = solver
    chunks = Vector{Any}(undef, P)
    U = Vector{typeof(u0)}(undef, P+1)
    T = Vector{eltype(tspan)}(undef, P+1)
    F = similar(U)
    return TimeParallelIterate(chunks, U, T, F)
end

# ---------------------------------------------------------------------------- #
#                                   Functions                                  #
# ---------------------------------------------------------------------------- #

# solution[k][n] ≡ solution.iterates[k].chunks[n]

"""
    length(iterate::TimeParallelIterate)

returns the number of chunks of `iterate`.
"""
Base.length(iterate::TimeParallelIterate) = length(iterate.chunks)

"""
    getindex(iterate::TimeParallelIterate, n::Integer)

returns the `n`-th chunk of a [`TimeParallelIterate`](@ref).
"""
Base.getindex(iterate::TimeParallelIterate, n::Integer) = iterate.chunks[n]

"""
    setindex!(iterate::TimeParallelIterate, value, n::Integer)

stores `value` into the `n`-th chunk of a [`TimeParallelIterate`](@ref).
"""
Base.setindex!(iterate::TimeParallelIterate, value, n::Integer) = iterate.chunks[n] = value

"""
    lastindex(iterate::TimeParallelIterate)

returns the last index of `iterate`.
"""
Base.lastindex(iterate::TimeParallelIterate) = lastindex(iterate.chunks)

"""
    show(io::IO, iterate::TimeParallelIterate)

prints a full description of `iterate` and its contents to a stream `io`.
"""
Base.show(io::IO, iterate::TimeParallelIterate) = NSDEBase._show(io, iterate)

"""
    summary(io::IO, iterate::TimeParallelIterate)

prints a brief description of `iterate` to a stream `io`.
"""
Base.summary(io::IO, iterate::TimeParallelIterate) = NSDEBase._summary(io, iterate)

# ---------------------------------------------------------------------------- #
#                                    Methods                                   #
# ---------------------------------------------------------------------------- #

function (iterate::TimeParallelIterate)(t::Real)
    N = length(iterate)
    if t ≤ iterate[1].t[1]
        return iterate[1].t[1], iterate[1].u[1]
    elseif t ≥ iterate[end].t[end]
        return iterate[end].t[end], iterate[end].u[end]
    else
        for n = 1:N
            if iterate[n].t[1] ≤ t ≤ iterate[n].t[end]
                return iterate[n](t)
            end
        end
    end
end
