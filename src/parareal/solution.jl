# src/parareal/solution.jl

"""
    PararealSolution <: AbstractTimeParallelSolution

A composite type for an [`AbstractTimeParallelSolution`](@ref) obtained using [`Parareal`](@ref).

# Constructors
```julia
PararealSolution(lastiterate, errors)
PararealSolution(problem::AbstractInitialValueProblem, parareal::Parareal)
```

## Arguments
- `lastiterate :: PararealIterate`
- `errors :: AbstractVector{ℝ} where ℝ<:Real` : iteration errors.

# Functions
- [`firstindex`](@ref) : first index.
- [`getindex`](@ref) : get iterate.
- [`lastindex`](@ref) : last index.
- [`numiterates`](@ref) : number of iterates.
- [`numchunks`](@ref) : number of chunks of last iterate.
- [`setindex!`](@ref) : set iterate.

# Methods

    (solution::PararealSolution)(t::Real)

returns the value of `solution` at `t` via interpolation.
"""
mutable struct PararealSolution{
            lastiterate_T<:PararealIterate,
            errors_T<:AbstractVector{<:Real},
            iterates_T<:Union{AbstractVector{<:PararealIterate}, Nothing},
        } <: AbstractTimeParallelSolution
    lastiterate::lastiterate_T
    errors::errors_T
    iterates::iterates_T
end

function PararealSolution(problem::AbstractInitialValueProblem, parareal::Parareal; saveiterates::Bool=false)
    lastiterate = PararealIterate(problem, parareal)
    @↓ K = parareal.parameters
    @↓ ϵ_T ← typeof(ϵ) = parareal.tolerance
    errors = Vector{ϵ_T}(undef, K)
    iterates = saveiterates ? [PararealIterate(problem, parareal) for i in 1:K] : nothing
    return PararealSolution(lastiterate, errors, iterates)
end

# ---------------------------------- METHODS ----------------------------------

(solution::PararealSolution)(tₚ::Real) = solution.lastiterate(tₚ)

# --------------------------------- FUNCTIONS ---------------------------------

"""
    length(solution::PararealSolution)

returns the number of chunks of `solution`.
"""
Base.length(solution::PararealSolution) = length(solution.lastiterate)

"""
    numchunks(solution::PararealSolution)

returns the number of chunks of `solution`.
"""
numchunks(solution::PararealSolution) = numchunks(solution.lastiterate)

"""
    numiterates(solution::PararealSolution)

returns the number of iterates of `solution`.
"""
numiterates(solution::PararealSolution) = length(solution.errors)

"""
    getindex(solution::PararealSolution, n::Integer)

returns the `n`-th chunk of the last iteration of a [`PararealSolution`](@ref).
"""
Base.getindex(solution::PararealSolution, n::Integer) = solution.lastiterate[n]

"""
    setindex!(solution::PararealSolution, chunk::AbstractInitialValueSolution, n::Integer)

stores an `AbstractInitialValueSolution` as the `n`-th chunk of the last iteration of a [`PararealSolution`](@ref).
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
    TimeParallelSolution(problem::AbstractInitialValueProblem, parareal::Parareal)

returns a [`PararealSolution`](@ref) constructor for the solution of `problem` with `parareal`.
"""
TimeParallelSolution(problem::AbstractInitialValueProblem, parareal::Parareal) = PararealSolution(problem, parareal)

"""
    collect_iterates!(solution::PararealSolution; directory)

reads the per-iterate chunk files written by the MPI backend under
`saveiterates` back into `solution.iterates`. Chunks below the diagonal
(`n < k`) are final from earlier iterations and are copied forward.
"""
function collect_iterates!(solution::PararealSolution; directory::String)
    solution.iterates isa Nothing && throw(ArgumentError("`collect_iterates!` needs a solution built with `saveiterates = true`."))
    for k = 1:numiterates(solution)
        for n = 1:k-1
            solution.iterates[k][n] = solution.iterates[k-1][n]
        end
        for n = k:numchunks(solution)
            filename = joinpath(directory, "iter_$(k)_chunk_$(n).jls")
            if isfile(filename)
                local_data = open(deserialize, filename, "r")
                if local_data.chunk_n !== nothing
                    solution.iterates[k][n] = local_data.chunk_n
                end
            end
        end
    end
    return solution
end

"""
    flatten(solution::PararealSolution)

[`flatten`](@ref)s the last iterate: one `(u, t)` pair for the whole span,
chunk boundaries deduplicated.
"""
flatten(solution::PararealSolution) = flatten(solution.lastiterate)

"""
    seams(solution::PararealSolution) :: Vector{<:Real}

[`seams`](@ref) of the last iterate.
"""
seams(solution::PararealSolution) = seams(solution.lastiterate)

"""
    maxseam(solution::PararealSolution) :: Real

[`maxseam`](@ref) of the last iterate.
"""
maxseam(solution::PararealSolution) = maxseam(solution.lastiterate)
