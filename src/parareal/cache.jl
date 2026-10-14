"""
    PararealCache <: AbstractTimeParallelCache

Pre-allocated state for one [`Parareal`](@ref) run: the chunk boundary values
`U`, the fine and coarse boundary results `F` and `G`, the previous-iterate
snapshot `U_` for the error function, the chunk time grid `T`, the chunk
problems (built once, their `u0` slots referencing `U` so the correction
updates them in place), and reusable solver caches. Method-dependent fields can
be replaced when switching solver, RHS, or weight-rule types; the state shape
and chunk count stay fixed.
Each run refreshes the initial state from the supplied problem. Moving a cache
to a different span requires `shiftwindow!`. RHS changes rebuild the chunk caches;
do not mutate RHS internals in place while retaining their solver caches.
"""
mutable struct PararealCache{
    makeGs_T<:AbstractVector{Bool},
    U_T<:AbstractVector{<:AbstractVector{<:Number}},
    T_T<:AbstractVector{<:Real},
} <: AbstractTimeParallelCache
    makeGs::makeGs_T
    U::U_T
    T::T_T
    F::U_T
    G::U_T
    U_::U_T
    chunkproblems
    finecaches
    coarsecache
    coarsesolution
    finesolver
    coarsesolver
    weightstate
    F_::U_T
    duration::Float64
    settings
    chunk_failed::Vector{Bool}
end

function PararealCache(
    problem::AbstractInitialValueProblem,
    parareal::Parareal;
    duration=(problem.tspan[2] - problem.tspan[1]) / parareal.parameters.N,
)
    @↓ u0, (t0, tN) ← tspan = problem
    @↓ finesolver, coarsesolver = parareal
    @↓ N = parareal.parameters
    validate_duration(t0, tN, N, duration)
    finesolver = chunk_solver(finesolver, duration)
    coarsesolver = chunk_solver(coarsesolver, duration)

    makeGs = trues(N)
    makeGs[1] = false

    T = [(N - n + 1) / N * t0 + (n - 1) / N * tN for n in 1:(N + 1)] # stable sum

    # U[1] is a COPY of the initial value, not an alias: window drivers (MoWi)
    # inject fresh starts into U[1] with copyto!, which must never write into
    # the user's u0.
    U = Vector{typeof(u0)}(undef, N)
    U[1] = copy(u0)
    for n in 2:N
        U[n] = similar(u0)
    end
    F = [n == 1 ? copy(u0) : similar(u0) for n in 1:N]
    G = [n == 1 ? copy(u0) : similar(u0) for n in 1:N]
    U_ = [similar(u0) for n in 1:N]

    # Chunk initial states alias the U buffers updated by the correction.
    chunkproblems = [remake(problem; u0=U[n], tspan=(T[n], T[n + 1])) for n in 1:N]

    finecaches = [NSDEBase.initialize_cache(chunkproblems[n], finesolver) for n in 1:N]
    coarsecache = NSDEBase.initialize_cache(chunkproblems[1], coarsesolver)
    coarsesolution = NSDEBase.initialize_solution(chunkproblems[1], coarsesolver)

    return PararealCache(
        makeGs,
        U,
        T,
        F,
        G,
        U_,
        chunkproblems,
        finecaches,
        coarsecache,
        coarsesolution,
        finesolver,
        coarsesolver,
        WeightState(weight_rule(parareal.criterion)),
        [similar(u0) for _ in 1:N],
        float(duration),
        (parareal.finesolver, parareal.coarsesolver, float(duration)),
        fill(false, N),
    )
end

#---------------------------------- FUNCTIONS ----------------------------------

function TimeParallelCache(problem::AbstractInitialValueProblem, parareal::Parareal; kwargs...)
    return PararealCache(problem, parareal; kwargs...)
end

"""
    shiftwindow!(cache::PararealCache, τ0, τN)

re-targets a cache at the window `[τ0, τN]`: recomputes the chunk grid `T` and
updates every chunk problem's `tspan` in place. This is the official seam for
window drivers (e.g. NSDEMovingWindow) — chunk problems are built once at
cache construction, so their time spans must be moved through this function,
never by poking `cache.T` alone.
"""
function shiftwindow!(cache::PararealCache, τ0::Real, τN::Real; duration=cache.duration)
    validate_duration(τ0, τN, length(cache.U), duration)
    cache.duration = duration
    @↓ T, chunkproblems = cache
    N = length(chunkproblems)
    for n in 1:(N + 1)
        T[n] = (N - n + 1) / N * τ0 + (n - 1) / N * τN # stable sum
    end
    for n in 1:N
        tspan = (T[n], T[n + 1])
        @↑ chunkproblems[n] = tspan
    end
    return cache
end

"""
    prepare_chunks!(cache::PararealCache, parareal::Parareal)

rebuilds the chunk schedules and solver caches after a change of window or method.
"""
function prepare_chunks!(cache::PararealCache, parareal::Parareal)
    parareal.parameters.N == length(cache.U) ||
        throw(DimensionMismatch("cache chunk count differs; build a new cache when changing N"))
    duration = cache.duration
    settings = (parareal.finesolver, parareal.coarsesolver, duration)
    settings === cache.settings && return cache
    cache.finesolver = chunk_solver(parareal.finesolver, duration)
    cache.coarsesolver = chunk_solver(parareal.coarsesolver, duration)
    cache.finecaches = [
        NSDEBase.initialize_cache(chunk, cache.finesolver) for chunk in cache.chunkproblems
    ]
    cache.coarsecache = NSDEBase.initialize_cache(cache.chunkproblems[1], cache.coarsesolver)
    cache.coarsesolution = NSDEBase.initialize_solution(cache.chunkproblems[1], cache.coarsesolver)
    cache.settings = settings
    return cache
end

"""
    reference(problem, parareal::Parareal) :: Vector

returns the serial fine solution at the `N + 1` chunk interfaces, computed on
exactly the same absolute-time chunk grid and schedule as `parareal`.
"""
function reference(problem::AbstractInitialValueProblem, parareal::Parareal)
    N = parareal.parameters.N
    t0, tN = problem.tspan
    T = [(N - n + 1) / N * t0 + (n - 1) / N * tN for n in 1:(N + 1)]
    fine = chunk_solver(parareal.finesolver, (tN - t0) / N)
    values = [copy(problem.u0)]
    for n in 1:N
        chunk = remake(problem; u0=copy(values[end]), tspan=(T[n], T[n + 1]))
        sol = NSDEBase.solve(chunk, fine)
        push!(values, copy(chunk_endpoint(sol, T[n + 1])))
    end
    return values
end

"""
    chunk_endpoint(sol, t)

returns the value of `sol` at the chunk endpoint `t`: the stored final value when
the solver stopped exactly at `t`, otherwise the interpolant (for generic solvers
that overshoot).
"""
chunk_endpoint(sol, t) = sol.t[end] == t ? sol.u[end] : sol(t)

# Construction and window shifts share the same nominal-duration contract.
function validate_duration(t0, tN, N, duration)
    isfinite(t0) && isfinite(tN) && tN > t0 ||
        throw(ArgumentError("time span must be finite and increasing"))
    isfinite(duration) && duration > 0 ||
        throw(ArgumentError("chunk duration must be finite and positive"))
    abs((tN - t0) - N * duration) <= 32eps(max(abs(float(t0)), abs(float(tN)), N * duration)) ||
        throw(ArgumentError("window span differs from N times the nominal chunk duration"))
    return nothing
end

"""
    prepare_problem!(cache::PararealCache, problem, parareal)

refreshes the initial data and RHS caches from `problem`; moving the time grid
needs [`shiftwindow!`](@ref).
"""
function prepare_problem!(cache::PararealCache, problem, parareal)
    parareal.parameters.N == length(cache.U) ||
        throw(DimensionMismatch("cache chunk count differs; build a new cache when changing N"))
    problem.tspan == (first(cache.T), last(cache.T)) || throw(
        ArgumentError(
            "problem time span differs from cache; call shiftwindow! or build a new cache"
        ),
    )
    axes(problem.u0) == axes(cache.U[1]) && eltype(problem.u0) == eltype(cache.U[1]) ||
        throw(DimensionMismatch("cache state shape or element type differs; build a new cache"))
    fill!(cache.chunk_failed, false)
    copyto!(cache.U[1], problem.u0)
    cache.makeGs[1] = false
    if problem.rhs !== cache.chunkproblems[1].rhs
        cache.chunkproblems = [
            remake(problem; u0=cache.U[n], tspan=(cache.T[n], cache.T[n + 1])) for
            n in eachindex(cache.U)
        ]
        cache.settings = nothing
    end
    return cache
end
