struct PararealCache{
            makeGs_T<:AbstractVector{<:Bool},
            U_T<:AbstractVector{<:AbstractVector{<:Number}},
            T_T<:AbstractVector{<:Real},
            F_T<:AbstractVector{<:AbstractVector{<:Number}},
            G_T<:AbstractVector{<:AbstractVector{<:Number}},
            U__T<:AbstractVector{<:AbstractVector{<:Number}},
        } <: AbstractTimeParallelCache
    makeGs::makeGs_T
    U::U_T
    T::T_T
    F::F_T
    G::G_T
    U_::U__T
end

function PararealCache(problem::AbstractInitialValueProblem, parareal::Parareal)
    @↓ u0, (t0, tN) ← tspan = problem
    @↓ N = parareal.parameters
    makeGs = trues(N)
    U = Vector{typeof(u0)}(undef, N)
    U[1] = u0
    makeGs[1] = false
    F = similar(U)
    G = similar(U)
    U_ = similar(U)
    T = Vector{typeof(t0)}(undef, N+1)
    for n = 1:N+1
        T[n] = (N - n + 1) / N * t0 + (n - 1) / N * tN # stable sum
    end
    return PararealCache(makeGs, U, T, F, G, U_)
end

#---------------------------------- FUNCTIONS ----------------------------------

TimeParallelCache(problem::AbstractInitialValueProblem, parareal::Parareal) = PararealCache(problem, parareal)
