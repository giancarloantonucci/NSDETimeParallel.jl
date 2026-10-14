# Guided tour: Parareal

Parareal partitions the interval into `N` chunks. A sequential coarse pass guesses
the chunk starts; a parallel batch of fine solves evaluates those guesses. A
coarse correction updates the starts, and another fine batch evaluates them.
Stopping early saves parallel rounds, but changes the accuracy achieved.

Run this page in order or use `examples/guided_tour.jl`. Start Julia with
`--threads=4` to exercise the threaded example.

## First run and its status

```@example parallel
using NSDERungeKutta, NSDETimeParallel, LinearAlgebra
problem = Logistic(0.3, (0.0, 2.0))
fine, coarse = RK4(h=0.005), RK4(h=0.05)
par = Parareal(fine, coarse; N=4,
    criterion=RelativeUpdate(atol=1e-12, rtol=1e-9))
r = solve(problem, par; history=:interfaces)
@assert r.diagnostics.status in (:converged, :exact)
@assert r.diagnostics.fine_batches == r.diagnostics.sweeps + 1
@assert r.errors === r.diagnostics.residuals
r.diagnostics
```

The default backend is serial. Use it to establish a baseline before parallel
execution. The default criterion is `RelativeUpdate()`: every stored chunk start
must pass a mixed absolute/relative update test. Its first test follows a
correction, so the first residual is `Inf`. That sentinel alone is not failure.
Always read `diagnostics.status`:

| Bulk status | Meaning |
|---|---|
| `:converged` | The chosen criterion passed |
| `:exact` | The matching discrete fine reference was reached by exact completion |
| `:budget_exhausted` | The fine-batch budget ended before acceptance |
| `:nonfinite` | A numerical failure or nonfinite propagation prevented a usable result |

An allocated but unused solution has `:not_started`. A numerical failure before
the first fine batch leaves empty residual and weight histories.

## Verify exact completion first

```@example parallel
full = remake(par; criterion=NoStop())
e = solve(problem, full; history=:interfaces)
Ustar = reference(problem, full)
@assert e.diagnostics.status == :exact
@assert e.diagnostics.fine_batches == 4
@assert e.diagnostics.sweeps == 3
@assert boundaryvalues(e.lastiterate) == Ustar
for (k, U) in enumerate(e.diagnostics.interfaces)
    @assert U[1:k] == Ustar[1:k]
end
nothing # hide
```

There are `N` stored starts and `N+1` returned boundary values; the extra value
is the last fine endpoint. At batch `k`, at least the first `k` starts are exact,
and the fine solves from them reproduce the first `k` chunks. At `N` batches,
all chunks are exact, using `N-1` corrections. `NoStop()` disables tolerance
acceptance, not the budget or failure handling: a smaller `K` can stop it early.

Exactness is relative to the **chunked discrete fine reference**, with matching
schedules, absolute times, arithmetic, and deterministic finite propagations.
It is not exactness relative to the ODE. A plain serial solve with a different
step schedule need not be bitwise identical.

## An end-to-end accuracy and cost study

Check the coarse solver over a chunk, then compare final interfaces with the
matching fine reference. Do not integrate the whole long problem merely to
inspect its first chunk.

```@example parallel
chaos = Lorenz([2.0, 3.0, -14.0], (0.0, 2.0))
F = RK4(h=0.002)
chunk = remake(chaos; tspan=(0.0, 0.5))
fine_end = solve(chunk, F).u[end]
for hG in (0.05, 0.02, 0.01)
    coarse_end = solve(chunk, RK4(h=hG)).u[end]
    println((hG=hG, chunk_error=norm(coarse_end-fine_end)))
end
p = Parareal(F, RK4(h=0.02); N=4, criterion=RelativeUpdate(rtol=1e-8))
a = solve(chaos, p; backend=ThreadsBackend())
ref = reference(chaos, p)
interface_error = maximum(norm.(boundaryvalues(a.lastiterate) .- ref))
ζ = costratio(chaos, p)
model = theoretical_speedup(a, ζ)
budget = iteration_budget(2.0, 4, ζ)
@assert isfinite(interface_error) && isfinite(model)
@assert 0 <= budget <= 4
println((status=a.diagnostics.status, batches=a.diagnostics.fine_batches,
         modelled_speedup=model, interface_error, max_seam=maxseam(a)))
nothing # hide
```

`costratio` measures coarse/fine time on the first chunk after warm-up. Its value
is machine-dependent and can be noisy on tiny examples. `theoretical_speedup`
then applies a model; it does not measure parallel elapsed time.
For `K` bulk batches, `s=K-1` and
`K_G=K*(N-1)-K*(K-1)/2`. The default model is
`N/(K + ζ*K_G)`, assuming enough workers, equal chunk costs, and no communication.
Thus its ideal ceiling is `N/K` even with free coarse work. See
[Representation and accounting](conventions.md) for other models and limitations.

`iteration_budget` answers how many batches that cost model permits, not how many
will converge. Zero means the requested speed-up is infeasible in the model.
`contractionrate(a)` fits the positive finite residual history; it can return
`NaN` if there are fewer than two usable entries. It is not a certified derivative
bound, particularly when an adaptive weight changes the residual's meaning.

## Change the accuracy target deliberately

```@example parallel
weighted = remake(p; criterion=WeightedDefect(ϵ=1e-8,
    weight=FixedWeight(exp(0.9))))  # an illustrative rate per unit time
wrun = solve(chaos, weighted)
secant = remake(p; criterion=WeightedDefect(ϵ=1e-8, weight=SecantWeight()))
srun = solve(chaos, secant)
@assert !isempty(srun.diagnostics.weights)
measured_weight = last(srun.diagnostics.weights)
weighted_error = Wnorm(srun, reference(chaos, secant))
@assert isfinite(weighted_error) && measured_weight >= 1
maximum_rule = WeightedDefect(ϵ=1e-8, weight=FixedWeight(1), p=Inf)
nothing # hide
```

`FixedWeight(w)` is a **per-unit-time base**. A rate `λ` gives
`w=exp(max(0,λ))`; a per-chunk base `wc` converts as `exp(log(wc)/ΔT)`.
The rate above illustrates the convention; it is not a bound inferred for this run.
A secant weight measures amplification between successive solves of the same
chunk. It can rise or fall and is not a uniform Lipschitz bound. Read its measured
value from diagnostics, not `SecantWeight().w`.

`WeightedDefect(p=1)` averages discounted defects over all `N` stored starts,
including the zero initial defect. It can accept the coarse starts on the first
fine batch. The returned chunks are still **fine solves** from those starts.
The maximum form `p=Inf` avoids dilution by averaging, but still discounts late
jumps. Inspect `seams`, `maxseam`, and unweighted interface error as appropriate.
`Wnorm` uses all returned interfaces except the fixed initial value, including
the terminal fine output. Neither stopping rule certifies invariant statistics.

## Histories and cache reuse

```@example parallel
saved = solve(problem, par; history=:iterates)
@assert saved.iterates !== nothing
cache = initialize_cache(problem, par)
reused = initialize_solution(problem, par)
new_problem = remake(problem; u0=[0.4])
par(cache, reused, new_problem)
@assert boundaryvalues(reused.lastiterate)[1] == [0.4]
larger = remake(par; N=8)
@assert larger.parameters.K == 8
nothing # hide
```

A one-off `solve(...; history=:iterates)` allocates the required history itself.
For an existing buffer, allocate it with that option before requesting full history.
`:interfaces` stores boundary vectors after batches; `:none` keeps only the final
iterate and ordinary diagnostics. Full history can dominate memory use.

Reuse requires compatible state shape/type and chunk count, and the same time
span. Rebuild buffers when those change. Reusing a configuration resets its running
weight. `remake(par; N=...)` lets `K` follow `N` when the previous budget equalled
the old `N`; an explicitly smaller budget stays smaller.

## Experiments to try next

1. Increase coarse `h` and compare batch count, seam size, and modelled cost.
2. Increase the Lorenz horizon with `N` fixed, checking accuracy as well as counts.
3. Compare relative updates, unweighted defects, fixed discounting, and secants.
4. Repeat a `NoStop()` reference check after changing the solver or chunk schedule.

A poor result can mean unstable coarse propagation, inaccurate phase or amplitude,
long chunks, a tight criterion, or an expensive coarse method. Batch growth need
not follow a simple exponential law. Weighting changes acceptance; it does not
repair an unstable solver or prove faster pointwise convergence. Very cheap chunks
may offer no wall-clock gain once scheduling and communication are included.
