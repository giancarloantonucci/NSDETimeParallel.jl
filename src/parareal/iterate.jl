"""
    PararealIterate <: AbstractTimeParallelIterate

A composite type for a single iterate in a [`PararealSolution`](@ref).

# Constructors
```julia
PararealIterate(chunks::AbstractVector{𝕊}) where 𝕊<:AbstractInitialValueSolution
PararealIterate(problem::AbstractInitialValueProblem, parareal::Parareal)
```

# Functions
- [`firstindex`](@ref) : first index.
- [`getindex`  ](@ref) : get chunk.
- [`lastindex` ](@ref) : last index.
- [`length`    ](@ref) : number of chunks.
- [`setindex!` ](@ref) : set chunk.

# Methods

    (iterate::PararealIterate)(t::Real)
    
returns the value of `iterate` at `t` via interpolation.
"""
struct PararealIterate{chunks_T<:(AbstractVector{𝕊} where 𝕊<:AbstractInitialValueSolution)} <: AbstractTimeParallelIterate
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

function (iterate::PararealIterate)(tₚ::Real)
    N = length(iterate)
    if tₚ < iterate[1].t[1]
        return iterate[1](tₚ)
    end
    for n = 1:N
        if (n > 1 ? iterate[n-1].t[end] : iterate[n].t[1]) ≤ tₚ < iterate[n].t[end]
            return iterate[n](tₚ)
        end
    end
    if tₚ ≥ iterate[N].t[end]
        return iterate[N](tₚ)
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

returns the `n`-th chunk of `iterate`.
"""
Base.getindex(iterate::PararealIterate, n::Integer) = iterate.chunks[n]

"""
    setindex!(iterate::PararealIterate, value::AbstractInitialValueSolution, n::Integer)

stores `value` into the `n`-th chunk of `iterate`.
"""
Base.setindex!(iterate::PararealIterate, value::AbstractInitialValueSolution, n::Integer) = iterate.chunks[n] = value

"""
    firstindex(iterate::PararealIterate)

returns the first index of `iterate`.
"""
Base.firstindex(iterate::PararealIterate) = firstindex(iterate.chunks)

"""
    lastindex(iterate::PararealIterate)

returns the last index of `iterate`.
"""
Base.lastindex(iterate::PararealIterate) = lastindex(iterate.chunks)
