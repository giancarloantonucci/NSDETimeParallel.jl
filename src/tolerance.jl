"standard error function."
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
        r += norm(U[n] - V[n]) / norm(U[n])
    end
    U_previous .= U
    return r / N
end

"weighted error function."
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
