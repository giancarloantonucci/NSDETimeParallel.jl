# src/parareal/cache.jl

"""
    PararealCache <: AbstractTimeParallelCache

Pre-allocated state for one [`Parareal`](@ref) run: the chunk boundary values
`U`, the fine and coarse boundary results `F` and `G`, the previous-iterate
snapshot `U_` for the error function, the chunk time grid `T`, the chunk
problems (built once, their `u0` slots referencing `U` so the correction
updates them in place), and reusable solver caches.
"""
struct PararealCache{
            makeGs_T<:AbstractVector{Bool},
            U_T<:AbstractVector{<:AbstractVector{<:Number}},
            T_T<:AbstractVector{<:Real},
            chunkproblems_T<:AbstractVector{<:AbstractInitialValueProblem},
            finecaches_T<:AbstractVector{<:AbstractInitialValueCache},
            coarsecache_T<:AbstractInitialValueCache,
            coarsesolution_T<:AbstractInitialValueSolution,
        } <: AbstractTimeParallelCache
    makeGs::makeGs_T
    U::U_T
    T::T_T
    F::U_T
    G::U_T
    U_::U_T
    chunkproblems::chunkproblems_T
    finecaches::finecaches_T
    coarsecache::coarsecache_T
    coarsesolution::coarsesolution_T
end

function PararealCache(problem::AbstractInitialValueProblem, parareal::Parareal)
    @↓ u0, (t0, tN) ← tspan = problem
    @↓ finesolver, coarsesolver = parareal
    @↓ N = parareal.parameters

    makeGs = trues(N)
    makeGs[1] = false

    T = [(N - n + 1) / N * t0 + (n - 1) / N * tN for n = 1:N+1] # stable sum

    # U[1] is a COPY of the initial value, not an alias: window drivers (MoWi)
    # inject fresh starts into U[1] with copyto!, which must never write into
    # the user's u0.
    U = Vector{typeof(u0)}(undef, N)
    U[1] = copy(u0)
    for n = 2:N
        U[n] = similar(u0)
    end
    F = [n == 1 ? copy(u0) : similar(u0) for n = 1:N]
    G = [n == 1 ? copy(u0) : similar(u0) for n = 1:N]
    U_ = [similar(u0) for n = 1:N]

    # Chunk problems are built ONCE; their `u0` fields reference the U buffers,
    # which the correction updates with `copyto!` — so these stay valid for the
    # whole run and the per-chunk-per-iteration `copy(problem, …)` is gone.
    chunkproblems = [copy(problem, U[n], T[n], T[n+1]) for n = 1:N]

    finecaches = [NSDEBase.initialize_cache(chunkproblems[n], finesolver) for n = 1:N]
    coarsecache = NSDEBase.initialize_cache(chunkproblems[1], coarsesolver)
    coarsesolution = NSDEBase.initialize_solution(chunkproblems[1], coarsesolver)

    return PararealCache(makeGs, U, T, F, G, U_, chunkproblems, finecaches, coarsecache, coarsesolution)
end

#---------------------------------- FUNCTIONS ----------------------------------

TimeParallelCache(problem::AbstractInitialValueProblem, parareal::Parareal) = PararealCache(problem, parareal)

"""
    shiftwindow!(cache::PararealCache, τ0, τN)

re-targets a cache at the window `[τ0, τN]`: recomputes the chunk grid `T` and
updates every chunk problem's `tspan` in place. This is the official seam for
window drivers (e.g. NSDEMovingWindow) — chunk problems are built once at
cache construction, so their time spans must be moved through this function,
never by poking `cache.T` alone.
"""
function shiftwindow!(cache::PararealCache, τ0::Real, τN::Real)
    @↓ T, chunkproblems = cache
    N = length(chunkproblems)
    for n = 1:N+1
        T[n] = (N - n + 1) / N * τ0 + (n - 1) / N * τN # stable sum
    end
    for n = 1:N
        tspan = (T[n], T[n+1])
        @↑ chunkproblems[n] = tspan
    end
    return cache
end
