"Immutable weight-rule settings; bases are per unit time."
abstract type AbstractWeightRule end
"Use a fixed positive per-unit-time discount base."
struct FixedWeight{T<:Real} <: AbstractWeightRule
    w::T
    function FixedWeight(w::Real=1.0)
        isfinite(w) && w > 0 || throw(ArgumentError("weight must be finite and positive"))
        return new{typeof(float(w))}(float(w))
    end
end
"Latest-pair same-chunk secant rule with a scaled denominator guard and a floor of one."
struct SecantWeight{T<:Real} <: AbstractWeightRule
    w::T
    τ::T
    δ::T
end
function SecantWeight(; w::Real=1.0, τ::Real=1e-10, δ::Real=1.0)
    isfinite(w) && w > 0 || throw(ArgumentError("weight must be finite and positive"))
    isfinite(τ) && τ >= 0 || throw(ArgumentError("guard must be finite and nonnegative"))
    isfinite(δ) && δ > 0 || throw(ArgumentError("divisor must be finite and positive"))
    return SecantWeight(promote(float(w), float(τ), float(δ))...)
end
"Legacy neighbouring-chunk running maximum; prefer an explicit fixed or secant rule."
struct NeighbourWeight{T<:Real} <: AbstractWeightRule
    w::T
    δ::T
end
function NeighbourWeight(; w::Real=1.0, δ::Real=1.0)
    isfinite(w) && w > 0 || throw(ArgumentError("weight must be finite and positive"))
    isfinite(δ) && δ > 0 || throw(ArgumentError("divisor must be finite and positive"))
    return NeighbourWeight(promote(float(w), float(δ))...)
end
weight_rule(w::AbstractWeightRule) = w
mutable struct WeightState{R<:AbstractWeightRule,T<:Real}
    rule::R
    w::T
end
WeightState(rule::AbstractWeightRule) = WeightState(rule, rule.w)
needs_fine_history(::AbstractWeightRule) = false
needs_fine_history(::SecantWeight) = true
update_weight!(state::WeightState{<:FixedWeight}, cache, solution, k) = state
function update_weight!(state::WeightState{<:NeighbourWeight}, cache, solution, k)
    U, F, T = cache.U, cache.F, cache.T
    candidate = state.w
    for i in 2:(length(U) - 1)
        den = norm(U[i]-U[i - 1])
        den > eps(typeof(den))*max(norm(U[i]), norm(U[i - 1])) || continue
        ratio = norm(F[i + 1]-F[i])/den
        measured = ratio^(1/(T[i + 1]-T[i])) / state.rule.δ
        isfinite(measured) && (candidate=max(candidate, measured))
    end
    state.w = candidate
    return state
end
function update_weight!(state::WeightState{<:SecantWeight}, cache, solution, k)
    rule = state.rule
    if k > 1
        best = -Inf
        for n in eachindex(cache.U)
            u, v = cache.U[n], cache.U_[n]
            den = norm(u - v)
            den > rule.τ * max(1, norm(u), norm(v)) || continue
            ratio = norm(solution.lastiterate[n].u[end] - cache.F_[n]) / den
            isfinite(ratio) || continue
            candidate = ratio^(1 / (cache.T[n + 1] - cache.T[n])) / rule.δ
            isfinite(candidate) && (best = max(best, candidate))
        end
        isfinite(best) && (state.w = max(one(state.w), best))
    end
    for n in eachindex(cache.F_)
        copyto!(cache.F_[n], solution.lastiterate[n].u[end])
    end
    return state
end
