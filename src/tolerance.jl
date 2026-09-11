# src/tolerance.jl

"""
    ψ₁(cache, k, weights)

standard relative error function: the mean over chunk boundaries of
`‖U[n] − V[n]‖ / ‖U[n]‖`, where `V` is the previous iterate (the coarse
prediction `F` at `k = 1`). Ignores `weights`.
"""
function ψ₁(cache, k, weights)
    @↓ U, T, U_previous ← U_ = cache
    if k > 1
        V = U_previous
    else
        @↓ V ← F = cache
    end
    r = 0.0
    N = length(U)
    for n = 1:N
        d = norm(U[n])
        δ = norm(U[n] - V[n])
        r += iszero(d) ? δ : δ / d
    end
    for n = 1:N
        copyto!(U_previous[n], U[n]) # deep snapshot: `.=` would alias the inner arrays
    end
    return r / N
end

"""
    ψ₂(cache, k, weights)

weighted error function: the mean over chunk boundaries of
`‖w^(T[1] − T[n]) (U[n] − F[n])‖`. With `w = exp(Λ)` for a problem with
Lyapunov exponent `Λ`, this DISCOUNTS each boundary error at the rate the
dynamics amplifies it — the moving-window criterion of the thesis. On long
chaotic spans it deliberately forgives late-window errors, so convergence in
ψ₂ means "early boundaries settled", not "uniformly within ϵ"; check the
terminal error against a serial fine solve when that distinction matters
(e.g. in speed-up benchmarks, where K sets the Amdahl ceiling).
"""
function ψ₂(cache, k, weights)
    @↓ U, T, F = cache
    @↓ w = weights
    w = max(1.0, w)
    r = 0.0
    N = length(U)
    for n = 1:N
        Wₙ = w^(T[1] - T[n])
        r += norm(Wₙ * (U[n] - F[n]))
    end
    return r / N
end

"""
    ψ∞(cache, k, weights)

local (pipelined) proximity function — the ℓ∞ member of the weighted family
(thesis §3.4, `subsec:local_proximity`): the MAXIMUM over chunk boundaries of
`‖w^(T[1] − T[n]) (U[n] − F[n])‖`, against `ψ₂`'s mean. Same defect, same
per-unit-time base `w`, so `ψ₂ ≤ ψ∞ ≤ N ψ₂` on any state.

Why it exists — two properties `ψ₂` cannot have:

- **Sharper certification.** `ψ∞(U) ≤ ϵ` certifies the weighted-max error
  `max_n w^(T[1]−T[n]) ‖U[n] − U*[n]‖ ≤ s(θ_F) ϵ` with
  `s(θ_F) = Σ_{j<N} θ_F^j`, `θ_F = Λ_F/w` — no factor `N`
  (thesis `thm:equivalence_error_norm_local`; the ℓ¹ constant is `N·s(θ_F)`).
  To certify a target radius `r`, set `ϵ = r / s(θ_F)`.
- **Prefix-decomposable acceptance.** `ψ∞ ≤ ϵ` is the conjunction of
  per-chunk tests, decidable left to right with one boolean riding the
  messages a pipeline already sends. [`PipelinedMPIBackend`](@ref) therefore
  HONOURS `ψ∞` — the one criterion it can enforce without re-serialising the
  pipeline — while `ψ₁`/`ψ₂` remain bulk-synchronous only.

Like `ψ₂`, evaluated over the stored chunk STARTS (the terminal boundary is
omitted — see the convention footnote in the docs) and clamping `w ≥ 1`. In
the bulk-synchronous loop the `n = 1` term is identically zero
(`F[1]` mirrors `U[1]`).
"""
function ψ∞(cache, k, weights)
    @↓ U, T, F = cache
    @↓ w = weights
    w = max(1.0, w)
    r = 0.0
    N = length(U)
    for n = 1:N
        Wₙ = w^(T[1] - T[n])
        r = max(r, norm(Wₙ * (U[n] - F[n])))
    end
    return r
end

"""
    Tolerance <: AbstractTolerance

A composite type for the tolerance mechanism of an time-parallel solver.

# Constructors
```julia
Tolerance(ϵ, ψ, weights)
Tolerance(; ϵ=1e-12, ψ=ψ₁, weights=Weights())
```

## Arguments
- `ϵ :: Real` : tolerance.
- `ψ :: Function` : error function.
- `weights :: Weights` : weights for ψ.
"""
struct Tolerance{ϵ_T<:Real, ψ_T<:Function, weights_T<:Weights} <: AbstractTolerance
    ϵ::ϵ_T
    ψ::ψ_T
    weights::weights_T
end

Tolerance(; ϵ::Real=1e-12, ψ::Function=ψ₁, weights::Weights=Weights()) = Tolerance(ϵ, ψ, weights)
