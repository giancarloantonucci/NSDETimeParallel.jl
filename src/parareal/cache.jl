"""
    PararealCache <: AbstractTimeParallelCache

A composite type for the [`AbstractTimeParallelCache`](@ref) of [`Parareal`](@ref).

# Constructors
```julia
PararealCache(U, F, G, T)
PararealCache(problem::AbstractInitialValueProblem, parareal::Parareal)
```

# Arguments
- `U :: AbstractVector{<:Union{Number, AbstractVector{<:Number}}}` : chunks' initial conditions.
- `F :: AbstractVector{<:Union{Number, AbstractVector{<:Number}}}` : chunks' fine values.
- `G :: AbstractVector{<:Union{Number, AbstractVector{<:Number}}}` : chunks' coarse values.
- `T :: AbstractVector{<:Number}` : chunks' initial time.
"""
struct PararealCache{U_T, F_T, G_T, T_T} <: AbstractTimeParallelCache
    U::U_T
    F::F_T
    G::G_T
    T::T_T
end

function PararealCache(problem::AbstractInitialValueProblem, parareal::Parareal)
    @↓ P = parareal
    @↓ u0, (t0, tN) ← tspan = problem
    u0_T = typeof(u0)
    U = Vector{u0_T}(undef, P+1)
    U[1] = u0
    F = similar(U)
    G = similar(U)
    t0_T = typeof(t0)
    T = Vector{t0_T}(undef, P+1)
    T[1] = t0
    return PararealCache(U, F, G, T)
end
