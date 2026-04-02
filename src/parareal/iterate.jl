# src/parareal/iterate.jl

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
- [`getindex`](@ref) : get chunk.
- [`lastindex`](@ref) : last index.
- [`length`](@ref) : number of chunks.
- [`numchunks`](@ref) : number of chunks.
- [`setindex!`](@ref) : set chunk.

# Methods

    (iterate::PararealIterate)(t::Real)

returns the value of `iterate` at `t` via interpolation.
"""
struct PararealIterate{
            # U_T<:AbstractVector{<:AbstractVector{<:Number}},
            # T_T<:AbstractVector{<:Real},
            chunks_T<:(AbstractVector{𝕊} where 𝕊<:AbstractInitialValueSolution),
        } <: AbstractTimeParallelIterate
    # U::U_T
    # T::T_T
    chunks::chunks_T
end

# function PararealIterate(problem::AbstractInitialValueProblem, parareal::Parareal)
#     # @↓ u0, t0 ← tspan[1] = problem
#     @↓ N = parareal.parameters
#     # U = Vector{typeof(u0)}(undef, N)
#     # U[1] = u0
#     # T = Vector{typeof(t0)}(undef, N+1)
#     # T[1] = t0
#     chunks = Vector{AbstractInitialValueSolution}(undef, N)
#     return PararealIterate(chunks)
# end

function PararealIterate(problem::AbstractInitialValueProblem, parareal::Parareal)
    @↓ N = parareal.parameters
    chunks = Vector{AbstractInitialValueSolution}(undef, N)
    return PararealIterate(chunks)
end

#----------------------------------- METHODS -----------------------------------

function (iterate::PararealIterate)(tₚ::Real)
    N = length(iterate)
    # if tₚ < iterate[1].t[1]
    #     return iterate[1](tₚ)
    # end
    for n = 1:N
        if (n > 1 ? iterate[n-1].t[end] : iterate[n].t[1]) ≤ tₚ < iterate[n].t[end]
            return iterate[n](tₚ)
        end
    end
    if tₚ == iterate[N].t[end]
        return iterate[N](tₚ)
    end
    # if tₚ ≥ iterate[N].t[end]
    #     return iterate[N](tₚ)
    # end
end
# TODO: maybe add a flag iscollected to use one method or the other
# function (iterate::PararealIterate)(tₚ::Real)
#     @↓ U, T = iterate
#     idx = 0
#     for n = 1:length(T)-1
#         if T[n] ≤ tₚ ≤ T[n+1] # TODO: Maybe < T[n+1] to avoid unnecessary error? Check
#             idx = n
#             break
#         end
#     end
#     uₚ = @fetchfrom workers()[idx] NSDETimeParallel.chunkfinesolution(tₚ)
#     return uₚ
# end

#---------------------------------- FUNCTIONS ----------------------------------

"""
    length(iterate::PararealIterate)

returns the number of chunks of `iterate`.
"""
Base.length(iterate::PararealIterate) = length(iterate.chunks)

"""
    numchunks(iterate::PararealIterate)

returns the number of chunks of `iterate`.
"""
numchunks(iterate::PararealIterate) = length(iterate.chunks)

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

function Wnorm(iterate::PararealIterate, reference::AbstractInitialValueSolution, w::Number)
    N = length(iterate)
    Ts = [[iterate[n].t[begin] for n = 1:N]; iterate[N].t[end]]
    Wₙ(n) = w ^ (Ts[n] / (Ts[n-1] - Ts[n]))
    Uₙᵏ(n) = iterate(Ts[n])
    Uₙ⁺(n) = reference(Ts[n])
    return norm([Wₙ(n) * (Uₙᵏ(n) - Uₙ⁺(n)) for n = 2:N+1])
end
