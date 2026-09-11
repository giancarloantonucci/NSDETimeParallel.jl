# src/utils.jl

"""
    theoretical_speedup(K::Integer, N::Integer, ζ::Real; wallclock=false, estimate=:work) :: Real

Parareal's own speed-up ceiling over the serial fine solve, for `K` iterations
over `N` chunks with cost ratio `ζ` = (coarse solve time over one chunk) /
(fine solve time over one chunk). Measure ζ — [`costratio`](@ref) does it on
one chunk of each — rather than guessing it. Measured speed-ups must be judged
against this ceiling, not against the serial fine solve wishfully.

Three estimates (thesis §2.7), selected with `estimate`:

- `:work` (default) — the classical WORK bound, whose `(1 - (K - 1)/2N)`
  factor credits the shrinking `k:N` sweeps (`eq:improved_parallel_efficiency`).
  That credit is saved work, not saved wall time: with workers ≥ chunks a
  parallel sweep costs one chunk-time however few chunks remain, so wall-clock
  measurements at small `N` (and any `K = N` point) sit BELOW this bound by
  construction.
- `:wallclock` — the work bound with the credit dropped, `N / (K (1 + ζN))`:
  the right comparator for a timing benchmark with a full worker pool.
  `wallclock = true` is an equivalent spelling, kept for backward
  compatibility; an explicit `estimate` wins if both are given.
- `:aubanel` — the task-scheduled (pipelined) ceiling
  (`eq:aubanel_parallel_efficiency`),
  `N / (K (1 - (K - 1)/2N)(1 + ζ) + ζ(N - 1))`: the serial coarse chain is
  passed rank to rank and hidden behind fine work, so the `ζN` term no longer
  multiplies `K`. This is the comparator for [`PipelinedMPIBackend`](@ref).
  At `ζ = 0` it coincides with `:work` — with a free coarse solver the
  pipeline has nothing to hide.
"""
function theoretical_speedup(K::Integer, N::Integer, ζ::Real; wallclock::Bool=false,
                             estimate::Symbol=(wallclock ? :wallclock : :work))
    K ≥ 1 && N ≥ 1 && ζ ≥ 0 || throw(ArgumentError("`theoretical_speedup` needs K ≥ 1, N ≥ 1, ζ ≥ 0."))
    credit = 1 - (K - 1) / (2N)
    estimate === :work      && return N / (K * (1 + ζ * N) * credit)
    estimate === :wallclock && return N / (K * (1 + ζ * N))
    estimate === :aubanel   && return N / (K * credit * (1 + ζ) + ζ * (N - 1))
    throw(ArgumentError("`theoretical_speedup` estimate must be :work, :wallclock or :aubanel, got :$estimate."))
end

"""
    iteration_budget(S::Real, N::Integer, ζ::Real; estimate=:improved) :: Int

the largest per-window iteration count `K ∈ 1:N` for which Parareal over `N`
chunks with cost ratio `ζ` still meets the target speed-up `S`, from the
parallel-efficiency estimates behind [`theoretical_speedup`](@ref); `0` when
the target is infeasible — no `K ≥ 1` reaches `S`, not even a single sweep.
The budget is capped at `N` because the algorithm terminates after `N`
sweeps: a target so modest that every sweep count meets it yields `N`, not a
number beyond the algorithm's range.

- `estimate = :improved` (default): the work-credited bound,
  `S(K) = N / (K (1 + ζN)(1 − (K − 1)/2N))`, which is strictly decreasing on
  `1:N`. It is inverted by integer bisection against
  [`theoretical_speedup`](@ref) itself, so the budget agrees with the
  estimate it is derived from bit for bit. (The closed-form root of the
  underlying quadratic was used before; its discriminant goes negative
  exactly when the target is BELOW the `K = N` speed-up — every `K` feasible —
  and the old code read that as "infeasible" and returned `0`. It also
  cancelled badly at large `N`.)
- `estimate = :naive`: `K = ⌊N/(S(1 + ζN))⌋`, likewise clamped to `0:N`.

Measure `ζ` with [`costratio`](@ref); a step-size guess distorts the budget
exactly where it matters (large `ζN`).
"""
function iteration_budget(S::Real, N::Integer, ζ::Real; estimate::Symbol=:improved)
    S > 0 && N ≥ 1 && ζ ≥ 0 || throw(ArgumentError("`iteration_budget` needs S > 0, N ≥ 1, ζ ≥ 0."))
    if estimate === :naive
        # Clamp in floating point BEFORE converting: a tiny S makes the quotient
        # far larger than typemax(Int), and `floor(Int, ·)` would throw.
        return Int(clamp(floor(N / (S * (1 + ζ * N))), 0.0, float(N)))
    elseif estimate === :improved
        speedup(K) = theoretical_speedup(K, N, ζ; estimate=:work)
        speedup(1) < S && return 0 # not even one sweep reaches the target
        speedup(N) ≥ S && return N # every admissible sweep count does
        # speedup is decreasing on 1:N, speedup(lo) ≥ S > speedup(hi): bisect.
        lo, hi = 1, N
        while hi - lo > 1
            mid = (lo + hi) ÷ 2
            if speedup(mid) ≥ S
                lo = mid
            else
                hi = mid
            end
        end
        return lo
    else
        throw(ArgumentError("`iteration_budget` estimate must be :improved or :naive, got :$estimate."))
    end
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
iteration exhibits no observable rate, and downstream bounds (e.g. the §4.3
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
function costratio(problem::AbstractInitialValueProblem, parareal::AbstractTimeParallelSolver; repeats::Integer=3)
    repeats ≥ 1 || throw(ArgumentError("`costratio` needs repeats ≥ 1."))
    @↓ finesolver, coarsesolver = parareal
    @↓ u0, (t0, tN) ← tspan = problem
    N = parareal.parameters.N
    chunk = copy(problem, u0, t0, t0 + (tN - t0) / N)
    NSDEBase.solve(chunk, finesolver)   # compile warm-up
    NSDEBase.solve(chunk, coarsesolver)
    tF = tG = Inf
    for _ = 1:repeats
        tF = min(tF, @elapsed NSDEBase.solve(chunk, finesolver))
        tG = min(tG, @elapsed NSDEBase.solve(chunk, coarsesolver))
    end
    return tG / tF
end
