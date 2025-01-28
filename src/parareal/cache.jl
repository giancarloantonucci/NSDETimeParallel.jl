struct PararealCache{
            skips_T<:AbstractVector{<:Bool},
            F_T<:AbstractVector{<:AbstractVector{<:Number}},
            G_T<:AbstractVector{<:AbstractVector{<:Number}},
            U_T<:AbstractVector{<:AbstractVector{<:Number}},
        } <: AbstractTimeParallelCache
    skips::skips_T
    F::F_T
    G::G_T
    U_::U_T
end

function PararealCache(problem::AbstractInitialValueProblem, parareal::Parareal)
    @↓ N = parareal.parameters
    skips = falses(N)
    skips[1] = true
    @↓ u0, t0 ← tspan[1] = problem
    F = Vector{typeof(u0)}(undef, N)
    G = similar(F)
    U_ = similar(F)
    return PararealCache(skips, F, G, U_)
end

#---------------------------------- FUNCTIONS ----------------------------------

TimeParallelCache(problem::AbstractInitialValueProblem, parareal::Parareal) = PararealCache(problem, parareal)
