"""
    Weights <: AbstractWeights

A composite type for the weights of [`Tolerance`](@ref).

# Constructors
```julia
Weights(; w=1.0, updatew=false)
```

## Arguments
- `w :: Union{AbstractVector{ℝ}, ℝ} where ℝ<:Real` : weighting factor for ψ.
- `updatew :: Bool` : flags when to [`update!`](@ref) `w` using (an approximation of) the Lipschitz function of the fine solver.

# Functions
- [`update!`](@ref) : updates `w` using (an approximation of) the Lipschitz function of the fine solver.
"""
mutable struct Weights{w_T<:(Union{AbstractVector{ℝ}, ℝ} where ℝ<:Real), updatew_T<:Bool} <: AbstractWeights
    w::w_T
    updatew::updatew_T
end

Weights(; w::Union{AbstractVector{ℝ}, ℝ}=1.0, updatew::Bool=false) where ℝ<:Real = Weights(w, updatew)

#---------------------------------- FUNCTIONS ----------------------------------

"""
    update!(weights::Weights, U<:AbstractVector{𝕍}, F<:AbstractVector{𝕍}) where 𝕍<:AbstractVector{ℂ} where ℂ<:Number

updates `weights.w` based on `U` and `F`.
"""
function update!(weights::Weights, U::AbstractVector{𝕍}, F::AbstractVector{𝕍}) where 𝕍<:AbstractVector{ℂ} where ℂ<:Number
    @↓ w, updatew = weights
    # TODO: Add `a` in Weights for Adaptive MoWi
    # @↓ w, updatew, a = weights
    a = 1
    N = length(U)
    w₁ = 0.0
    w₂ = 0.0
    if updatew
        for i = 2:N-1
            r = a * norm(F[i+1] - F[i]) / norm(U[i] - U[i-1])
            w₁ += r
            w₂ = max(w₂, r)
        end
        w₁ /= N - 1
    end
    a₁ = 0.0
    a₂ = 1.0
    w = max(w, a₁ * w₁ + a₂ * w₂)
    @↑ weights = w
    return weights
end
