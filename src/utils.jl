# src/utils.jl

function check_cost_inputs(N::Integer, ζ::Real)
    N >= 1 && isfinite(ζ) && ζ >= 0 ||
        throw(ArgumentError("cost model needs N ≥ 1 and a finite ζ ≥ 0"))
    return nothing
end

function speedup_from_counts(K_F, K_G, fine_calls, N, ζ, estimate)
    fine_cost = if estimate === :wallclock
        K_F
    elseif estimate === :work
        fine_calls / N
    else
        throw(
            ArgumentError(
                "estimate must be :wallclock or :work; pipeline speed-up requires its task schedule",
            ),
        )
    end
    return N / (fine_cost + ζ * K_G)
end

"""
    theoretical_speedup(K::Integer, N::Integer, ζ::Real; estimate=:wallclock) :: Real
    theoretical_speedup(diagnostics::PararealDiagnostics, ζ::Real; estimate=:wallclock) :: Real
    theoretical_speedup(solution::PararealSolution, ζ::Real; estimate=:wallclock) :: Real

Modelled speed-up over N serial fine chunks. `K` is the number of FINE BATCHES,
not corrections: a finite bulk return has `s = K-1`, with `1 ≤ K ≤ N`.
The cost ratio `ζ` is coarse/fine time per chunk; see [`costratio`](@ref).

The exact-prefix implementation uses
`K_G = K*(N-1) - K*(K-1)/2` coarse calls and
`K*N - K*(K-1)/2` fine calls, including the initial coarse pass and final fine
batch. The diagnostics overload uses the recorded call counts directly.

- `:wallclock` (default) charges one fine-chunk time per batch and sequential
  coarse calls: `N/(K + ζ*K_G)`. It assumes at least N workers and equal chunk
  costs, and excludes communication, synchronization, and load imbalance.
- `:work` divides the total fine-call cost by N, retaining the same sequential
  coarse cost. This assumes perfect redistribution of the shrinking fine work;
  it is not the critical path of the bulk backend.

These are cost models, not measured speed-ups or unconditional performance
bounds. The compatibility keyword `wallclock=false` selects `:work`, while an
explicit `estimate` takes precedence. Pipeline diagnostics, nonfinite or unstarted
runs, and runs with no fine batch are refused. A pipeline needs its task schedule;
the old `:aubanel` formula is no longer offered as a model of native execution.
"""
function theoretical_speedup(
    K::Integer,
    N::Integer,
    ζ::Real;
    wallclock::Bool=true,
    estimate::Symbol=(wallclock ? :wallclock : :work),
)
    check_cost_inputs(N, ζ)
    1 <= K <= N || throw(ArgumentError("K must count fine batches in 1:N"))
    # Convert before multiplication, avoiding integer overflow in the counts.
    k, n = float(K), float(N)
    triangular = k * (k - 1) / 2
    K_G = k * (n - 1) - triangular
    fine_calls = k * n - triangular
    return speedup_from_counts(K, K_G, fine_calls, N, ζ, estimate)
end

function theoretical_speedup(
    diagnostics::PararealDiagnostics,
    ζ::Real;
    wallclock::Bool=true,
    estimate::Symbol=(wallclock ? :wallclock : :work),
)
    d = diagnostics
    d.fine_batches === nothing && throw(
        ArgumentError(
            "pipeline diagnostics have no global fine batches; a task schedule is required"
        ),
    )
    d.status in (:converged, :exact, :budget_exhausted) &&
    d.fine_batches >= 1 &&
    d.initial_starts !== nothing || throw(
        ArgumentError("speed-up needs a finite completed bulk run with at least one fine batch")
    )
    N = length(d.initial_starts)
    check_cost_inputs(N, ζ)
    1 <= d.fine_batches <= N && d.coarse_chunk_calls >= 0 && d.fine_chunk_calls >= 1 ||
        throw(ArgumentError("invalid diagnostic call counts"))
    return speedup_from_counts(
        d.fine_batches, d.coarse_chunk_calls, d.fine_chunk_calls, N, ζ, estimate
    )
end

function theoretical_speedup(solution::PararealSolution, ζ::Real; kwargs...)
    return theoretical_speedup(solution.diagnostics, ζ; kwargs...)
end

"""
    iteration_budget(S::Real, N::Integer, ζ::Real; estimate=:wallclock) :: Int

Largest FINE-BATCH budget `K ∈ 1:N` that reaches target modelled speed-up S;
returns zero if none does. A returned K permits K-1 corrections. This is a cost
budget, not a prediction that the iteration converges within that budget.

Integer bisection inverts [`theoretical_speedup`](@ref) with the SAME cost model
and exact-prefix coarse count. Both `:wallclock` (default) and `:work` are
supported. The older names `:naive` and `:improved` remain aliases for these two
models, respectively; they now use the corrected counts too.
"""
function iteration_budget(S::Real, N::Integer, ζ::Real; estimate::Symbol=:wallclock)
    check_cost_inputs(N, ζ)
    isfinite(S) && S > 0 || throw(ArgumentError("target speed-up must be finite and positive"))
    model = if estimate === :naive
        :wallclock
    elseif estimate === :improved
        :work
    else
        estimate
    end
    model in (:wallclock, :work) || throw(ArgumentError("estimate must be :wallclock or :work"))
    speedup(K) = theoretical_speedup(K, N, ζ; estimate=model)
    speedup(1) < S && return 0
    speedup(N) >= S && return N
    lo, hi = 1, N
    while hi - lo > 1
        mid = lo + (hi - lo) ÷ 2
        if speedup(mid) >= S
            lo = mid
        else
            hi = mid
        end
    end
    return lo
end

"""
    contractionrate(errors::AbstractVector{<:Real}) :: Float64
    contractionrate(solution::AbstractTimeParallelSolution) :: Float64

the observed per-iteration contraction rate `β` of an error trace: the
least-squares slope of `log(errors[k])` against `k`, exponentiated, so that
`errors[k] ≈ C βᵏ`. Non-finite entries and EXACT ZEROS are excluded before
fitting — finite termination drives the last error to bitwise zero, which is
a property of the algorithm's fixed point, not of its rate. Returns `NaN`
when fewer than two usable points remain: a run that converges in one
iteration exhibits no observable rate, and downstream bounds (e.g. the
outer radius `R = r/βᴷ`) must treat that as unbounded rather than fabricate a
number.
"""
function contractionrate(errors::AbstractVector{<:Real})
    ks = [k for k in eachindex(errors) if isfinite(errors[k]) && errors[k] > 0]
    length(ks) < 2 && return NaN
    x = float.(ks)
    y = log.(Float64.(errors[ks]))
    x̄ = sum(x) / length(x)
    ȳ = sum(y) / length(y)
    slope = sum((x .- x̄) .* (y .- ȳ)) / sum(abs2, x .- x̄)
    return exp(slope)
end

contractionrate(solution::AbstractTimeParallelSolution) = contractionrate(solution.errors)

"""
    costratio(problem, parareal; repeats=3) :: Real

measures ζ, the coarse/fine cost ratio over ONE chunk of `problem` under
`parareal`'s chunking — the input [`theoretical_speedup`](@ref) needs. Both
solvers are timed on the first chunk with `@elapsed` (best of `repeats`,
after a compile warm-up). Measure ζ rather than guessing it from step-size
ratios: those ignore per-step cost differences between the two solvers.
"""
function costratio(
    problem::AbstractInitialValueProblem, parareal::AbstractTimeParallelSolver; repeats::Integer=3
)
    repeats ≥ 1 || throw(ArgumentError("`costratio` needs repeats ≥ 1."))
    @↓ finesolver, coarsesolver = parareal
    @↓ u0, (t0, tN) ← tspan = problem
    N = parareal.parameters.N
    chunk = remake(problem; u0, tspan=(t0, t0 + (tN - t0) / N))
    finesolver = chunk_solver(finesolver, (tN-t0)/N)
    coarsesolver = chunk_solver(coarsesolver, (tN-t0)/N)
    fcache = initialize_cache(chunk, finesolver)
    gcache = initialize_cache(chunk, coarsesolver)
    fsolution = initialize_solution(chunk, finesolver)
    gsolution = initialize_solution(chunk, coarsesolver)
    NSDEBase.solve!(fcache, fsolution, chunk, finesolver)   # compile warm-up
    NSDEBase.solve!(gcache, gsolution, chunk, coarsesolver)
    tF = tG = Inf
    for _ in 1:repeats
        tF = min(tF, @elapsed NSDEBase.solve!(fcache, fsolution, chunk, finesolver))
        tG = min(tG, @elapsed NSDEBase.solve!(gcache, gsolution, chunk, coarsesolver))
    end
    return tG / tF
end
