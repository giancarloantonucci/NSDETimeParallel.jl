struct PararealCache{skips_T<:AbstractVector{<:Bool}, U_T<:AbstractVector{<:AbstractVector{<:Number}}, F_T<:AbstractVector{<:AbstractVector{<:Number}}, G_T<:AbstractVector{<:AbstractVector{<:Number}}, T_T<:AbstractVector{<:Real}} <: AbstractTimeParallelCache
    skips::skips_T
    U::U_T
    F::F_T
    G::G_T
    T::T_T
end

function PararealCache(problem::AbstractInitialValueProblem, parareal::Parareal)
    @↓ N = parareal.parameters
    @↓ u0, t0 ← tspan[1] = problem
    skips = falses(N)
    skips[1] = true
    U = Vector{typeof(u0)}(undef, N)
    F = similar(U)
    G = similar(U)
    F[1] = G[1] = U[1] = u0
    T = Vector{typeof(t0)}(undef, N+1)
    T[1] = t0
    return PararealCache(skips, U, F, G, T)
end

#---------------------------------- FUNCTIONS ----------------------------------

TimeParallelCache(problem::AbstractInitialValueProblem, parareal::Parareal) = PararealCache(problem, parareal)
