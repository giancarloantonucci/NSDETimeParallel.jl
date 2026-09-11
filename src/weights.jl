# src/weights.jl

"""
    Weights <: AbstractWeights

A composite type for the weights of [`Tolerance`](@ref).

# Constructors
```julia
Weights(; w=1.0, updatew=false, δ=1.0)
```

## Arguments
- `w :: Real` : weighting factor for ψ, a PER-UNIT-TIME base (`w = exp(Λ)` for per-time rate Λ). A finite, positive scalar. (Earlier versions admitted a vector here; nothing downstream — `ψ₂`, `ψ∞`, [`update!`](@ref), the MoWi zoom — ever defined what a per-chunk vector of bases meant, and each threw a `MethodError` on one. A vector is now refused at construction.)
- `updatew :: Bool` : flags when to [`update!`](@ref) `w` using (an approximation of) the Lipschitz function of the fine solver. Bulk-synchronous backends only: under [`PipelinedMPIBackend`](@ref) with `ψ∞` there is no root that sees all boundaries at a common sweep, so the updater cannot run — that path REJECTS `updatew = true`; preset `w` (a known rate, or a frozen probe measurement — see the criterion docs).
- `δ :: Real` : safety divisor applied to the MEASURED rate only (`δ < 1`
  inflates the measured base — the thesis's safety factor `w = C·Λ̂` with
  `C = 1/δ`). It never touches the running `w`, so with `updatew = false`
  [`update!`](@ref) is a strict no-op regardless of `δ`.

# Functions
- [`update!`](@ref) : updates `w` using (an approximation of) the Lipschitz function of the fine solver.
"""
mutable struct Weights{w_T<:Real, updatew_T<:Bool, δ_T<:Real} <: AbstractWeights
    w::w_T
    updatew::updatew_T
    δ::δ_T
    # The invariants live in the inner constructor so that no entry point —
    # keyword, positional, or a future one — can store a NaN weight or a zero
    # divisor.
    function Weights(w::w_T, updatew::updatew_T, δ::δ_T) where {w_T<:Real, updatew_T<:Bool, δ_T<:Real}
        isfinite(δ) && δ > 0 || throw(ArgumentError("`Weights` needs a finite, positive safety divisor δ; got $δ."))
        isfinite(w) && w > 0 || throw(ArgumentError("`Weights` needs a finite, positive scalar w; got $w."))
        # Stored as a float: `update!` and the MoWi zoom assign measured or
        # rescaled values back into this field, which an integer type would
        # reject with an InexactError.
        wf = float(w)
        return new{typeof(wf), updatew_T, δ_T}(wf, updatew, δ)
    end
end

# A vector is refused with a message that says why, rather than a MethodError
# from the first ψ that touches it.
Weights(w::AbstractVector, updatew, δ) = throw(ArgumentError("`Weights` takes a single per-unit-time base `w::Real`; a per-chunk vector of bases has no defined meaning in ψ₂, ψ∞ or `update!`. Got a $(length(w))-vector."))

Weights(; w::Union{Real,AbstractVector}=1.0, updatew::Bool=false, δ::Real=1.0) = Weights(w, updatew, δ)

#---------------------------------- FUNCTIONS ----------------------------------

"""
    update!(weights::Weights, U, F, T)

updates `weights.w` from the chunk-boundary values — bulk-synchronous
backends only (the pipelined `ψ∞` path has no global view and rejects
`updatew = true`): the largest measured
per-chunk amplification `‖F[i+1] − F[i]‖ / ‖U[i] − U[i−1]‖`, converted to a
PER-UNIT-TIME base via `r^(1/(T[i+1] − T[i]))`. The conversion is not
optional: [`ψ₂`](@ref) consumes `w` as `w^(T[1] − T[n])` — an exponent in
TIME — so `w` must be `exp(Λ)` with `Λ` a per-time rate. Storing the raw
per-chunk ratio (the old behaviour) is right only when chunks are exactly one
time unit long, as in the Lorenz thesis set-up; with 10-unit chunks the
discount runs at ten times the honest rate, and ψ₂ goes blind past the first
boundary.
"""
function update!(weights::Weights, U::AbstractVector{𝕍}, F::AbstractVector{𝕍}, T::AbstractVector{<:Real}) where 𝕍<:AbstractVector{ℂ} where ℂ<:Number
    @↓ w, updatew, δ = weights
    updatew || return weights # strict no-op when frozen
    # TODO: Add `a` in Weights for Adaptive MoWi
    # @↓ w, updatew, a = weights
    N = length(U)
    w₁ = 0.0
    w₂ = 0.0
    m = 0 # number of usable ratios
    for i = 2:N-1
        den = norm(U[i] - U[i-1])
        num = norm(F[i+1] - F[i])
        # An unresolved ratio carries no information about the rate: two
        # coincident starts (an equilibrium, an already-exact prefix, a
        # converged run) give 0/0 = NaN, and `max(w, NaN)` is NaN — which then
        # blinds ψ₂ for the rest of the run. Skip those, and skip any ratio
        # that does not come out finite.
        scale = max(norm(U[i]), norm(U[i-1]))
        den > eps(typeof(den)) * scale || continue
        r = (num / den)^(1 / (T[i+1] - T[i])) # per-chunk amplification → per-time base
        isfinite(r) || continue
        w₁ += r
        w₂ = max(w₂, r)
        m += 1
    end
    m == 0 && return weights # nothing measurable this sweep: keep the previous weight
    w₁ /= m
    a₁ = 0.0
    a₂ = 1.0
    # δ divides the MEASURED candidate only. Dividing the running w (the old
    # behaviour, `max(w, ·)/δ`) turned any δ ≠ 1 into a per-iteration ratchet:
    # `update!` runs once per Parareal sweep, so a frozen weight decayed
    # geometrically and Zoom's adjustments were silently rescaled every k.
    candidate = (a₁ * w₁ + a₂ * w₂) / δ
    isfinite(candidate) || return weights
    w = max(w, candidate)
    @↑ weights = w
    return weights
end
