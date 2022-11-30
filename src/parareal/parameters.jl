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
end

PararealParameters(; N::Integer=10, K::Integer=N) = PararealParameters(N, K)
