"""
    PararealSolution <: AbstractTimeParallelSolution

A composite type for an [`AbstractTimeParallelSolution`](@ref) obtained using [`Parareal`](@ref).

# Constructors
```julia
PararealSolution(errors, lastiterate, iterates)
PararealSolution(problem::AbstractInitialValueProblem, parareal::Parareal)
```

## Arguments
- `errors :: AbstractVector{ℝ} where ℝ<:Real` : iteration errors.
- `lastiterate :: PararealIterate`
- `iterates :: AbstractVector{𝕊} where 𝕊<:PararealIterate`

# Functions
- [`firstindex`](@ref) : first index.
- [`getindex`](@ref) : get iterate.
- [`lasterror`](@ref) : last iteration error.
- [`lastindex`](@ref) : last index.
- [`length`](@ref) : number of iterates.
- [`setindex!`](@ref) : set iterate.

# Methods

    (solution::PararealSolution)(t::Real)

returns the value of `solution` at `t` via interpolation.
"""
mutable struct PararealSolution{errors_T<:(AbstractVector{ℝ} where ℝ<:Real), lastiterate_T<:PararealIterate, iterates_T<:Union{AbstractVector{𝕊} where 𝕊<:PararealIterate, Nothing}} <: AbstractTimeParallelSolution
    errors::errors_T
    lastiterate::lastiterate_T
    iterates::iterates_T
end

function PararealSolution(problem::AbstractInitialValueProblem, parareal::Parareal)
    @↓ saveiterates = parareal
    @↓ K = parareal.parameters
    @↓ ϵ_T ← typeof(ϵ) = parareal.tolerance
    errors = Vector{ϵ_T}(undef, K)
    lastiterate = PararealIterate(problem, parareal)
    iterates = saveiterates ? [PararealIterate(problem, parareal) for i in 1:K] : nothing
    return PararealSolution(errors, lastiterate, iterates)
end

#----------------------------------- METHODS -----------------------------------

(solution::PararealSolution)(tₚ::Real) = solution.lastiterate(tₚ)

#---------------------------------- FUNCTIONS ----------------------------------

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

"""
    ktherror(solution::PararealSolution, k::Integer)

returns the convergence error for the `k`-th iterate of `solution`.
"""
ktherror(solution::PararealSolution, k::Integer) = solution.errors[k]

"""
    lasterror(solution::PararealSolution)

returns the convergence error for the last iterate of `solution`.
"""
lasterror(solution::PararealSolution) = solution.errors[end]

"""
    TimeParallelSolution(problem::AbstractInitialValueProblem, parareal::Parareal)

returns a [`PararealSolution`](@ref) constructor for the solution of `problem` with `parareal`.
"""
TimeParallelSolution(problem::AbstractInitialValueProblem, parareal::Parareal) = PararealSolution(problem, parareal)

function Wnorm(solution::PararealSolution, reference::AbstractInitialValueSolution, w::Number)
    K = length(solution.errors)
    if solution.iterates isa Nothing
        return error("`Wnorm`: `(solution::PararealSolution).iterates` contains `nothing`. Solve with `saveiterates = true`.")
    else
        return [Wnorm(solution.iterates[k], reference, w) for k = 1:K]
    end
end
