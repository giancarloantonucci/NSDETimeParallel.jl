function 𝜑₁(U, F, T)
    r = 0.0
    for n = 2:length(U)
        r += norm(U[n] - F[n]) / norm(U[n])
    end
    return r
end

function 𝜑₂(U, F, T, λ = 0.0)
    r = 0.0
    for n = 2:length(U)
        Wₙ = exp(λ * (T[1] - T[n]))
        r += norm(Wₙ * (U[n] - F[n]))^2
    end
    return r
end

"""
    ErrorCheck

A composite type for the error control mechanism used by a [`TimeParallelSolver`](@ref).

# Constructors
```julia
ErrorCheck(; 𝜑 = 𝜑₁, ϵ = 1e-12)
```

# Arguments
- `𝜑 :: Function` : error control function.
- `ϵ :: Real` : tolerance

# Functions
- [`show`](@ref) : shows name and contents.
- [`summary`](@ref) : shows name.
"""
struct ErrorCheck{𝜑_T, ϵ_T}
    𝜑::𝜑_T
    ϵ::ϵ_T
end

ErrorCheck(; 𝜑 = 𝜑₁, ϵ = 1e-12) = ErrorCheck(𝜑, ϵ)

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
