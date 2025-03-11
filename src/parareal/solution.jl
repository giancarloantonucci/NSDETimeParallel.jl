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
    TimeParallelSolution(problem::AbstractInitialValueProblem, parareal::Parareal)

returns a [`PararealSolution`](@ref) constructor for the solution of `problem` with `parareal`.
"""
TimeParallelSolution(problem::AbstractInitialValueProblem, parareal::Parareal) = PararealSolution(problem, parareal)

# function Wnorm(solution::PararealSolution, reference::AbstractInitialValueSolution, w::Number)
#     K = length(solution.errors)
#     if solution.iterates isa Nothing
#         return error("`Wnorm`: `(solution::PararealSolution).iterates` contains `nothing`. Solve with `saveiterates = true`.")
#     else
#         return [Wnorm(solution.iterates[k], reference, w) for k = 1:K]
#     end
# end

# function collect!(solution::PararealSolution)
#     N = length(solution.lastiterate)
#     for n = 1:N
#         solution[n] = @fetchfrom workers()[n] NSDETimeParallel.chunkfinesolution
#     end
#     return solution
# end

function collect!(solution::PararealSolution; directory::String="results")
    for n = 1:numchunks(solution)
        filename = joinpath(directory, "lastiter_chunk_$(n).jls")
        if isfile(filename)
            open(filename, "r") do file
                local_data = deserialize(file)
                if local_data.chunk_n !== nothing
                    solution.lastiterate[n] = local_data.chunk_n
                end
            end
        end
    end
    return solution
end

function collect_iterates(iterates::AbstractVector{<:PararealIterate}; directory::String="results")
    @↓ iterates = solution
    for k = 1:numiterates(solution)
        for n = 1:numchunks(solution)
            filename = joinpath(directory, "iter_$(k)_chunk_$(n).jls")
            if isfile(filename)
                open(filename, "r") do file
                    local_data = deserialize(file)
                    if local_data.chunk_n !== nothing
                        iterates[k][n] = local_data.chunk_n
                    end
                end
            end
        end
    end
    return iterates
end

# function collect_iterates!(solution::PararealSolution; dir::String="results")
#     for k = 1:numiterates(solution)
#         for n = 1:k-1
#             solution.iterates[k][n] = solution.iterates[k-1][n]
#         end
#         for n = k:numchunks(solution)
#             filename = joinpath(dir, "iter_$(k)_chunk_$(n).jls")
#             if isfile(filename)
#                 open(filename, "r") do file
#                     local_data = deserialize(file)
#                     if local_data.chunk_n !== nothing
#                         solution.iterates[k][n] = local_data.chunk_n
#                     end
#                 end
#             end
#         end
#     end
#     return solution
# end
