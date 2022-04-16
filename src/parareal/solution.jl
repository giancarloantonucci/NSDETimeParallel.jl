"""
    PararealSolution <: AbstractTimeParallelSolution

A composite type for an [`AbstractTimeParallelSolution`](@ref) obtained using [`Parareal`](@ref).

# Constructors
```julia
PararealSolution(errors, lastiterate, alliterates[, saveiterates])
PararealSolution(problem::AbstractInitialValueProblem, parareal::Parareal; saveiterates::Bool=false)
```

## Arguments
- `errors :: AbstractVector{ℝ} where ℝ<:Real` : iteration errors.
- `lastiterate :: PararealIterate`
- `alliterates :: AbstractVector{𝕊} where 𝕊<:PararealIterate`
- `saveiterates :: Bool` : flag to fill `alliterates`.

# Functions
- [`firstindex`](@ref) : first index.
- [`getindex`](@ref) : get iterate.
- [`lastindex`](@ref) : last index.
- [`length`](@ref) : number of iterates.
- [`setindex!`](@ref) : set iterate.

# Methods

    (solution::PararealSolution)(t::Real)
    
returns the value of `solution` at `t` via interpolation.
"""
mutable struct PararealSolution{errors_T<:(AbstractVector{ℝ} where ℝ<:Real), lastiterate_T<:PararealIterate, alliterates_T<:Union{AbstractVector{𝕊} where 𝕊<:PararealIterate, Nothing}, saveiterates_T<:Bool} <: AbstractTimeParallelSolution
    errors::errors_T
    lastiterate::lastiterate_T
    alliterates::alliterates_T
    saveiterates::saveiterates_T
end

function PararealSolution(problem::AbstractInitialValueProblem, parareal::Parareal; saveiterates::Bool=false)
    @↓ K = parareal
    @↓ ϵ_T ← typeof(ϵ) = parareal.tolerance
    errors = Vector{ϵ_T}(undef, K)
    lastiterate = PararealIterate(problem, parareal)
    alliterates = saveiterates ? [PararealIterate(problem, parareal) for i in 1:K] : nothing
    return PararealSolution(errors, lastiterate, alliterates, saveiterates)
end

#####
##### Methods
#####

(solution::PararealSolution)(tₚ::Real) = solution.lastiterate(tₚ)

#####
##### Functions
#####

"""
    length(solution::PararealSolution)

returns the number of chunks of `solution`.
"""
Base.length(solution::PararealSolution) = length(solution.lastiterate)

"""
    getindex(solution::PararealSolution, n::Integer)

returns the `n`-th chunk of the last iteration of a [`PararealSolution`](@ref).
"""
Base.getindex(solution::PararealSolution, n::Integer) = solution.lastiterate[n]

"""
    setindex!(solution::PararealSolution, chunk::AbstractInitialValueSolution, n::Integer)

stores an [`AbstractInitialValueSolution`](@ref) as the `n`-th chunk of the last iteration of a [`PararealSolution`](@ref).
"""
Base.setindex!(solution::PararealSolution, chunk::AbstractInitialValueSolution, n::Integer) = solution.lastiterate[n] = chunk

"""
    firstindex(solution::PararealSolution)

returns the first index of `solution`.
"""
Base.firstindex(solution::PararealSolution) = firstindex(solution.lastiterate)

"""
    lastindex(solution::PararealSolution)

returns the last index of `solution`.
"""
Base.lastindex(solution::PararealSolution) = lastindex(solution.lastiterate)

TimeParallelSolution(problem::AbstractInitialValueProblem, parareal::Parareal; saveiterates::Bool=false) = PararealSolution(problem, parareal; saveiterates=saveiterates)
