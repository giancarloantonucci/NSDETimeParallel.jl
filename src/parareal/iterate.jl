"""
    PararealIterate <: AbstractTimeParallelIterate

A composite type for the single iterations within a [`PararealSolution`](@ref).

# Constructors
```julia
PararealIterate(chunks)
PararealIterate(problem::AbstractInitialValueProblem, parareal::Parareal)
```

# Arguments
- `chunks :: AbstractVector{<:AbstractInitialValueSolution}` : vector of chunk solutions.

# Functions
- [`getindex`](@ref) : get chunk.
- [`lastindex`](@ref) : last index.
- [`length`](@ref) : number of chunks.
- [`setindex!`](@ref) : set chunk.
- [`show`](@ref) : shows name and contents.
- [`summary`](@ref) : shows name.
"""
struct PararealIterate{chunks_T} <: AbstractTimeParallelIterate
    chunks::chunks_T
end

function PararealIterate(problem::AbstractInitialValueProblem, parareal::Parareal)
    @↓ P = parareal
    chunks = Vector{AbstractInitialValueSolution}(undef, P)
    return PararealIterate(chunks)
end

#####
##### Methods
#####

function (iterate::PararealIterate)(t::Real)
    N = length(iterate)
    T0 = iterate[1].t[1]
    TN = iterate[end].t[end]
    if t ≤ T0 || t ≥ TN
        return error("t = $t is out of bounds $((T0, TN)).")
    else
        for n = 1:N
            t0 = iterate[n].t[1]
            tN = iterate[n].t[end]
            if t0 ≤ t ≤ tN
                return iterate[n](t)
            end
        end
    end
end

#####
##### Functions
#####

# solution[k][n] ≡ solution.iterates[k].chunks[n]

"""
    length(iterate::PararealIterate)

returns the number of chunks of `iterate`.
"""
Base.length(iterate::PararealIterate) = length(iterate.chunks)

"""
    getindex(iterate::PararealIterate, n::Integer)

returns the `n`-th chunk of a [`PararealIterate`](@ref).
"""
Base.getindex(iterate::PararealIterate, n::Integer) = iterate.chunks[n]

"""
    setindex!(iterate::PararealIterate, value, n::Integer)

stores `value` into the `n`-th chunk of a [`PararealIterate`](@ref).
"""
Base.setindex!(iterate::PararealIterate, value, n::Integer) = iterate.chunks[n] = value

"""
    lastindex(iterate::PararealIterate)

returns the last index of `iterate`.
"""
Base.lastindex(iterate::PararealIterate) = lastindex(iterate.chunks)

# function RungeKutta.extract(iterate::PararealIterate, i::Integer)
#     N = length(iterate)
#     @↓ u1 ← u, t1 ← t = RungeKutta.extract(iterate[1], i)
#     u = eltype(u1)[]
#     t = eltype(t1)[]
#     append!(u, u1[1:end])
#     append!(t, t1[1:end])
#     for n = 2:N
#         @↓ u1 ← u, t1 ← t = RungeKutta.extract(iterate[n], i)
#         append!(u, u1[2:end])
#         append!(t, t1[2:end])
#     end
#     return RungeKutta.RungeKuttaSolution(u, t)
# end
