"""
    PararealIterate <: AbstractTimeParallelIterate

A composite type for a single iterate in a [`PararealSolution`](@ref): one
fine chunk solution per time chunk. The chunk storage is CONCRETE — chunks are
built eagerly from the fine solver, so `eltype(iterate.chunks)` is the fine
solver's concrete solution type, not an abstract box. This matters: these
objects cross thread and process boundaries in the hot path.

# Constructors
```julia
PararealIterate(chunks::AbstractVector{𝕊}) where 𝕊<:AbstractInitialValueSolution
PararealIterate(problem::AbstractInitialValueProblem, parareal::Parareal)
```

# Functions
- [`firstindex`](@ref)/[`getindex`](@ref)/[`lastindex`](@ref)/[`setindex!`](@ref) : chunk access
- [`length`](@ref)/[`numchunks`](@ref) : number of chunks

# Methods

    (iterate::PararealIterate)(t::Real)

returns the value of `iterate` at `t` via interpolation of the owning chunk.
"""
struct PararealIterate{chunks_T<:(AbstractVector{𝕊} where 𝕊<:AbstractInitialValueSolution)} <:
       AbstractTimeParallelIterate
    chunks::chunks_T
end

function PararealIterate(problem::AbstractInitialValueProblem, parareal::Parareal)
    @↓ u0, (t0, tN) ← tspan = problem
    @↓ finesolver = parareal
    @↓ N = parareal.parameters
    # One representative chunk problem fixes shape and type; every chunk spans
    # the same length, so all chunk solutions are structurally identical.
    duration = chunkduration(t0, tN, N)
    chunkproblem = remake(problem; u0, tspan=(gridtime(t0), gridtime(t0) + duration))
    finesolver = chunk_solver(finesolver, duration)
    chunks = [NSDEBase.initialize_solution(chunkproblem, finesolver) for _ in 1:N]
    return PararealIterate(chunks)
end

#----------------------------------- METHODS -----------------------------------

function (iterate::PararealIterate)(tₚ::Real)
    N = length(iterate)
    for n in 1:N
        if (n > 1 ? iterate[n - 1].t[end] : iterate[n].t[1]) ≤ tₚ < iterate[n].t[end]
            return iterate[n](tₚ)
        end
    end
    if tₚ ≥ iterate[N].t[end]
        return iterate[N](tₚ) # clamp above (the chunk clamps internally)
    end
    return iterate[1](tₚ) # clamp below
end

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
function Base.setindex!(iterate::PararealIterate, value::AbstractInitialValueSolution, n::Integer)
    return iterate.chunks[n] = value
end

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

"""
    boundarytimes(iterate::PararealIterate)

returns the `N + 1` chunk-boundary times of `iterate`: each chunk's start,
plus the final chunk's end. These are the points Parareal iterates on — the
grid every convergence measure ([`Wnorm`](@ref), per-iterate error studies)
is taken over.
"""
function boundarytimes(iterate::PararealIterate)
    return [[iterate[n].t[begin] for n in 1:length(iterate)]; iterate[end].t[end]]
end

"""
    boundaryvalues(iterate::PararealIterate)

returns the `N + 1` chunk-boundary values of `iterate`: each chunk's first
state, plus the final chunk's last state — the `U` vector of the Parareal
iteration, read off the chunk solutions. The returned states ALIAS the chunk
storage; copy them before mutating.
"""
function boundaryvalues(iterate::PararealIterate)
    return [[iterate[n].u[begin] for n in 1:length(iterate)]; [iterate[end].u[end]]]
end

"""
    flatten(iterate::PararealIterate)

concatenates the chunk solutions into one pair of vectors `(u, t)`, dropping
each chunk's FIRST point after chunk 1: boundaries are stored twice (chunk
`n`'s last point and chunk `n+1`'s corrected start), and keeping both would
double-weight every boundary in any downstream statistic. The kept value is
the fine endpoint; the dropped corrected start agrees with it only to the
solve tolerance. `t` is monotone.
"""
function flatten(iterate::PararealIterate)
    u = copy(iterate[1].u)
    t = copy(iterate[1].t)
    for n in 2:length(iterate)
        append!(u, @view iterate[n].u[2:end])
        append!(t, @view iterate[n].t[2:end])
    end
    return u, t
end

"""
    seams(iterate::PararealIterate) :: Vector{<:Real}

returns the `N − 1` chunk-boundary seam sizes of `iterate`: at each interior
boundary, `‖(start of chunk n+1) − (end of chunk n)‖` — the distance between
the corrected starting value `U[n+1]` and the fine endpoint `F(U[n])` it
should agree with at convergence. This is exactly the discontinuity that
[`flatten`](@ref) hides by dropping the corrected start, so any statistical
claim made on flattened output should report [`maxseam`](@ref) alongside it:
under the weighted criterion `ψ₂`, late-window seams are unconstrained BY
DESIGN (the discount forgives them), and only the seam sizes say which regime
a run actually exercised. Empty for a single-chunk iterate.
"""
function seams(iterate::PararealIterate)
    return [norm(iterate[n + 1].u[begin] - iterate[n].u[end]) for n in 1:(length(iterate) - 1)]
end

"""
    maxseam(iterate::PararealIterate) :: Real

returns the largest chunk-boundary seam of `iterate` (see [`seams`](@ref)),
or `0.0` for a single-chunk iterate.
"""
maxseam(iterate::PararealIterate) = length(iterate) > 1 ? maximum(seams(iterate)) : 0.0

"""
    Wnorm(iterate, reference, weight; base=:time)

Weighted Euclidean distance at the chunk interfaces. `reference` may be a
solution or the vector returned by `reference(problem, parareal)`. The base
is per unit time, matching all weight rules. `base=:chunk` explicitly selects
the per-chunk convention. Supply a scalar or `FixedWeight`; adaptive rule
objects contain only their initial settings and are refused. To use a measured
weight, pass `solution.diagnostics.weights[end]`, or call `Wnorm(solution, reference)`.
"""
function Wnorm(iterate::PararealIterate, reference, weight; base::Symbol=:time)
    base in (:time, :chunk) || throw(ArgumentError("base must be :time or :chunk"))
    weight isa Union{Real,FixedWeight} || throw(
        ArgumentError(
            "Wnorm needs a scalar or FixedWeight; use solution.diagnostics.weights[end] for an adaptive rule",
        ),
    )
    w = weight isa Real ? weight : weight.w
    isfinite(w) && w > 0 || throw(ArgumentError("weight must be finite and positive"))
    w = max(one(w), w)
    Ts = boundarytimes(iterate)
    ΔT = Ts[2] - Ts[1]
    U = boundaryvalues(iterate)
    reference isa AbstractVector &&
        length(reference) != length(U) &&
        throw(DimensionMismatch("reference must contain N+1 interfaces"))
    return norm([
        w^(-(Ts[n] - Ts[1]) / (base === :chunk ? ΔT : 1)) *
        (U[n] - (reference isa AbstractVector ? reference[n] : reference(Ts[n]))) for
        n in 2:length(U)
    ])
end

# Unsolved chunks must have both defined time spans and explicit missing values,
# including when a solution buffer still contains an earlier successful run.
function invalidate_chunks!(iterate::PararealIterate, T, indices)
    for n in indices
        invalidate_chunk!(iterate[n], (T[n], T[n + 1]))
    end
    return iterate
end

function invalidate_chunk!(chunk, tspan)
    resize!(chunk.t, 2)
    resize!(chunk.u, 2)
    chunk.t[1], chunk.t[2] = tspan
    chunk.u[1] = fill!(similar(chunk.u[1]), NaN)
    chunk.u[2] = copy(chunk.u[1])
    hasproperty(chunk, :k) && chunk.k !== nothing && empty!(chunk.k)
    return chunk
end
