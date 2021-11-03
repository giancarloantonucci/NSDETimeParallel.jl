function ψ₁(cache, solution, k, args...)
    @↓ U, T = cache
    # if k > 1
    #     V = getU(solution[k-1])
    #     # @↓ V ← U = getU(solution[k-1])
    # else
    @↓ V ← F = cache
    # end
    r = 0.0
    N = length(U)
    for n = 1:N
        r += norm(U[n] - V[n]) / norm(U[n])
    end
    return r / N
end

function ψ₂(cache, solution, k, weights)
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
    ErrorControl <: AbstractTimeParallelParameters

A composite type for the error control mechanism of an [`AbstractTimeParallelSolver`](@ref).

# Constructors
```julia
ErrorControl(; ϵ=1e-12, ψ=ψ₁, weights=ErrorWeights())
```

# Arguments
- `ϵ :: Real` : tolerance
- `ψ :: Function` : error function.
- `weights :: ErrorWeights` : weights for ψ.

# Functions
- [`show`](@ref) : shows name and contents.
- [`summary`](@ref) : shows name.
"""
struct ErrorControl{ϵ_T, ψ_T, weights_T} <: AbstractTimeParallelParameters
    ϵ::ϵ_T
    ψ::ψ_T
    weights::weights_T
end

ErrorControl(; ϵ=1e-12, ψ=ψ₁, weights=ErrorWeights()) = ErrorControl(ϵ, ψ, weights)
