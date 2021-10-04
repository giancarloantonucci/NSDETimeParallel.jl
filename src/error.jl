function 𝜑₁(solution, k, Λ)
    @↓ T = solution
    @↓ U = solution[k]
    k > 1 ? (@↓ W ← U = solution[k-1]) : (@↓ W ← F = solution[k])
    r = 0.0
    N = length(U)
    for n = 1:N
        r += norm(U[n] - W[n]) / norm(U[n])
    end
    return r / N
end

function 𝜑₂(solution, k, Λ)
    @↓ T = solution
    @↓ U, F = solution[k]
    Λ = max(1.0, Λ)
    r = 0.0
    N = length(U)
    for n = 1:N
        # Wₙ = exp(λ * (T[1] - T[n]))
        Wₙ = Λ ^ (T[1] - T[n])
        r += norm(Wₙ * (U[n] - F[n]))
    end
    return r / N
end

"""
    ErrorCheck

A composite type for the error control mechanism used by a [`TimeParallelSolver`](@ref).

# Constructors
```julia
ErrorCheck(; 𝜑 = 𝜑₁, ϵ = 1e-12, Λ = 1.0)
```

# Arguments
- `𝜑 :: Function` : error control function.
- `ϵ :: Real` : tolerance
- `Λ :: Real` : Lipschitz constant of fine solver.

# Functions
- [`show`](@ref) : shows name and contents.
- [`summary`](@ref) : shows name.
"""
struct ErrorCheck{𝜑_T, ϵ_T, Λ_T, updateΛ_T}
    𝜑::𝜑_T
    ϵ::ϵ_T
    Λ::Λ_T
    updateΛ::updateΛ_T
end

ErrorCheck(; 𝜑 = 𝜑₁, ϵ = 1e-12, Λ = 1.0, updateΛ = false) = ErrorCheck(𝜑, ϵ, Λ, updateΛ)

function update_Lipschitz(Λ, U, F)
    N = length(U)
    for i = 2:N-1
        tmp = norm(F[i+1] - F[i]) / norm(U[i] - U[i-1])
        Λ = max(Λ, tmp)
        # Λ += norm(F[i+1] - F[i]) / norm(U[i] - U[i-1])
    end
    # return max(1.0, Λ / (N-1))
    return max(1.0, Λ)
end

# ---------------------------------------------------------------------------- #
#                                   Functions                                  #
# ---------------------------------------------------------------------------- #

"""
    show(io::IO, error_check::ErrorCheck)

prints a full description of `error_check` and its contents to a stream `io`.
"""
Base.show(io::IO, error_check::ErrorCheck) = NSDEBase._show(io, error_check)

"""
    summary(io::IO, error_check::ErrorCheck)

prints a brief description of `error_check` to a stream `io`.
"""
Base.summary(io::IO, error_check::ErrorCheck) = NSDEBase._summary(io, error_check)
