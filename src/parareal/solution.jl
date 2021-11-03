"""
    PararealSolution <: AbstractTimeParallelSolution

A composite type for an [`AbstractTimeParallelSolution`](@ref) obtained using [`Parareal`](@ref).

# Constructors
```julia
PararealSolution(iterates, errors[, saveiterates])
PararealSolution(problem::AbstractInitialValueProblem, parareal::Parareal)
```

# Arguments
- `iterates :: AbstractVector{<:PararealIterate}` : initial value problem, e.g. an [``](@ref).
- `errors :: ` :
- `saveiterates :: ` :

# Functions
- [`getindex`](@ref) : get iterate.
- [`lastindex`](@ref) : last index.
- [`length`](@ref) : number of iterates.
- [`setindex!`](@ref) : set iterate.
- [`show`](@ref) : shows name and contents.
- [`summary`](@ref) : shows name.
"""
mutable struct PararealSolution{errors_T, lastiterate_T, alliterates_T, saveiterates_T} <: AbstractTimeParallelSolution
    errors::errors_T
    lastiterate::lastiterate_T
    alliterates::alliterates_T
    saveiterates::saveiterates_T
end

function PararealSolution(problem, parareal::Parareal; saveiterates::Bool=false)
    @↓ P, K = parareal
    @↓ ϵ = parareal.control
    ϵ_T = typeof(ϵ)
    errors = Vector{ϵ_T}(undef, K)
    lastiterate = saveiterates ? nothing : PararealIterate(problem, parareal)
    alliterates = saveiterates ? [PararealIterate(problem, parareal) for i in 1:K] : nothing
    return PararealSolution(errors, lastiterate, alliterates, saveiterates)
end

#####
##### Methods
#####

function (solution::PararealSolution)(t::Real)
    @↓ saveiterates = solution
    if saveiterates
        @↓ alliterates = solution
        return alliterates[end](t)
    else
        @↓ lastiterate = solution
        return lastiterate(t)
    end
end

function (solution::PararealSolution)(t::Real, n::Integer)
    @↓ saveiterates = solution
    if saveiterates
        @↓ alliterates = solution
        return alliterates[n](t)
    else
        println("WARNING: `saveiterates` is set to `false`, therefore no iteration except the last has been saved.")
        @↓ lastiterate = solution
        return lastiterate(t)
    end
end

# To-Do: create mask to select right chunk

#####
##### Functions
#####

"""
    length(solution::PararealSolution)

returns the number of (not necessarily saved) iterates of `solution`.
"""
Base.length(solution::PararealSolution) = length(solution.errors)

"""
    getindex(solution::PararealSolution, k::Integer)

returns the `k`-th iterate of a [`PararealSolution`](@ref), if saved; else, it returns the `k`-th chunk of the last iteration of a [`PararealSolution`](@ref).
"""
function Base.getindex(solution::PararealSolution, k::Integer)
    @↓ saveiterates = solution
    if saveiterates
        @↓ alliterates = solution
        return alliterates[k]
    else
        @↓ lastiterate = solution
        return lastiterate[k]
    end
end

"""
    setindex!(solution::PararealSolution, iterate::PararealIterate, k::Integer)

stores a [`PararealIterate`](@ref) as the `n`-th iterate of a [`PararealSolution`](@ref).
"""
function Base.setindex!(solution::PararealSolution, iterate::PararealIterate, k::Integer)
    return solution.alliterates[k] = iterate
end

"""
    setindex!(solution::PararealSolution, chunk::AbstractInitialValueSolution, n::Integer)

stores an [`AbstractInitialValueSolution`](@ref) as the `n`-th chunk of the last iteration of a [`PararealSolution`](@ref).
"""
function Base.setindex!(solution::PararealSolution, chunk::AbstractInitialValueSolution,
                        n::Integer)
    return solution.lastiterate[n] = chunk
end

"""
    lastindex(solution::PararealSolution)

returns the last index of `solution`.
"""
Base.lastindex(solution::PararealSolution) = lastindex(solution.alliterates)

# function RungeKutta.extract(solution::PararealSolution, i::Integer)
#     return RungeKutta.extract(solution[end], i)
# end

function Base.resize!(solution::PararealSolution, k::Integer)
    @↓ saveiterates = solution
    if saveiterates
        @↓ alliterates = solution
        resize!(alliterates, k)
        return solution
    end
end
