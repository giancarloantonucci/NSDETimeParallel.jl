"""
    ErrorWeights <: AbstractTimeParallelParameters

A composite type for the weights of [`ErrorControl`](@ref).

# Constructors
```julia
ErrorWeights(; w=1.0, updatew=false)
```

# Arguments
- `w::Union{Real, AbstractVector} = ` : weighting of ψ.
- `updatew :: Bool` : to update `w` using local information.

# Functions
- [`show`](@ref) : shows name and contents.
- [`summary`](@ref) : shows name.
"""
mutable struct ErrorWeights{w_T, updatew_T} <: AbstractTimeParallelParameters
    w::w_T
    updatew::updatew_T
end

ErrorWeights(; w=1.0, updatew=false) = ErrorWeights(w, updatew)

#####
##### Functions
#####

function update!(weights::ErrorWeights, U, F)
    @↓ w, updatew = weights
    if updatew
        N = length(U)
        for i = 2:N-1
            # w += norm(F[i+1] - F[i]) / norm(U[i] - U[i-1])
            w = max(w, norm(F[i+1] - F[i]) / norm(U[i] - U[i-1]))
        end
    end
    # w = max(1.0, w / (N-1))
    w = max(1.0, w)
    @↑ weights = w
    return weights
end
