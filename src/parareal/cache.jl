struct PararealCache{U_T<:(AbstractVector{𝕍} where 𝕍<:AbstractVector{ℂ} where ℂ<:Number), F_T<:(AbstractVector{𝕍} where 𝕍<:AbstractVector{ℂ} where ℂ<:Number), G_T<:(AbstractVector{𝕍} where 𝕍<:AbstractVector{ℂ} where ℂ<:Number), T_T<:(AbstractVector{ℂ} where ℂ<:Number)} <: AbstractTimeParallelCache
    U::U_T
    F::F_T
    G::G_T
    T::T_T
end
function PararealCache(problem::AbstractInitialValueProblem, parareal::Parareal)
    @↓ P = parareal
    @↓ u0, t0 ← tspan[1] = problem
    U = Vector{typeof(u0)}(undef, P)
    F = similar(U)
    G = similar(U)
    F[1] = G[1] = U[1] = u0
    T = Vector{typeof(t0)}(undef, P+1)
    T[1] = t0
    return PararealCache(U, F, G, T)
end

#####
##### Functions
#####

TimeParallelCache(problem::AbstractInitialValueProblem, parareal::Parareal) = PararealCache(problem, parareal)
