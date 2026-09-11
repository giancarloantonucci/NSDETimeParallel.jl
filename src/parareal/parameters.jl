# src/parareal/parameters.jl

"""
    PararealParameters <: AbstractPararealParameters

A composite type for the basic parameters of Parareal.

# Constructors
```julia
PararealParameters(N, K)
PararealParameters(; N=10, K=N)
```

## Arguments
- `N :: Integer` : number of time chunks/processors.
- `K :: Integer` : maximum number of iterations.
"""
mutable struct PararealParameters{N_T<:Integer, K_T<:Integer} <: AbstractPararealParameters
    N::N_T
    K::K_T
    function PararealParameters(N::N_T, K::K_T) where {N_T<:Integer, K_T<:Integer}
        N ≥ 1 && K ≥ 1 || throw(ArgumentError("PararealParameters needs N ≥ 1 and K ≥ 1, got N = $N, K = $K."))
        # Not an error — the run self-terminates — but every sweep past N is
        # provably wasted: finite termination makes iterate N exact.
        K > N && @warn "PararealParameters: K = $K exceeds N = $N; finite termination makes every sweep past N dead work."
        return new{N_T, K_T}(N, K)
    end
end

PararealParameters(; N::Integer=10, K::Integer=N) = PararealParameters(N, K)
