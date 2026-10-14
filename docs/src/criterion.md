# Stopping criteria

`Parareal(fine, coarse; N, criterion=...)` makes the acceptance target explicit.
A reported residual is a diagnostic of that criterion, not automatically an error
bound against the unknown serial fine trajectory.

## All-interface updates

`RelativeUpdate(atol=1e-12, rtol=1e-9)` is the default. Every stored start must satisfy

```math
\|U_n-V_n\| \le \mathrm{atol} + \mathrm{rtol}\max(\|U_n\|,\|V_n\|).
```

The residual is the maximum of the update divided by its tolerance. Acceptance
requires it to be at most one. The first test follows the first correction;
the initial batch has no previous corrected iterate. Small updates imply small
trajectory errors only with an appropriate convergence estimate.

`MeanRelativeUpdate(ϵ=...)` uses the mean relative update, including its absolute
fallback at zero norm. It does not require every interface to pass separately.
The iteration loop, not the criterion, owns the previous-iterate snapshot.

## Weighted defects

`WeightedDefect(ϵ=1e-9, weight=FixedWeight(1), p=1)` evaluates the normalized
p-mean of weighted interface defects. For `p=1` this is the averaged defect; for
`p=Inf` it is the maximum. States use the Euclidean norm. On stored starts,
`d[1]=0` and `d[n]=U[n]-F(U[n-1])` for `n=2:N`. The normalization includes all
`N` starts, so the averaged form divides by `N`, including the zero first defect.

Every rule uses a **per-unit-time base** `w`: the discount at `T[n]` is
`max(1,w)^(T[1]-T[n])`. A desired per-chunk base `wc` on chunks of length `ΔT`
is specified as `w=exp(log(wc)/ΔT)`. `Wnorm` uses the same convention. `Wnorm(solution, reference)` uses the last
measured weight. For an explicit weight, use `Wnorm(solution.lastiterate,
reference, w)` or `FixedWeight(w)`. Adaptive rule objects are refused here:
`SecantWeight().w` is an initial setting, not a measured weight. A run that fails
before its first finite fine batch has no measured weight.

Acceptance can occur on the initial coarse starts after the first fine batch.
The returned objects are fine trajectories from those starts. The terminal
output is supplied by the last fine solve and introduces no independent terminal
join. A separately corrected terminal unknown would require its own defect.

The mean can hide one large defect: `ψ₂≤ϵ` permits
`norm(d[n])≤N*ϵ*w^(T[n]-T[1])`. Increasing `N` or the weight can therefore make
acceptance easier without improving the tail. `p=Inf` removes the factor `N`,
but still discounts late defects. Inspect `seams` and `maxseam` when unweighted
continuity matters. Neither choice establishes statistical fidelity.

## Weight rules

- `FixedWeight(w)` keeps the configured base. For a prescribed Lyapunov rate
  `λ`, use `FixedWeight(exp(max(0,λ)))`.
- `SecantWeight(w=1, τ=1e-10, δ=1)` compares each chunk's input and fine output
  across the latest two batches. A denominator is usable when it exceeds
  `τ*max(1,norm(Unew),norm(Uold))`. The largest finite ratio, converted to a
  per-time base and divided by `δ`, sets the new base with floor one. It can
  fall. If no pair passes the guard, the previous base is retained.
- `NeighbourWeight(w=1, δ=1)` is the legacy neighbouring-chunk estimator. It
  retains a running maximum, so it can only increase the base.

Secants are observed amplifications, not uniform Lipschitz bounds. Both adaptive
rules keep their current value in the run cache; reusing a solver starts again
from its configured `w`. Retrieve a measured value from `solution.diagnostics.weights`.
For a later fixed-weight probe, construct a new criterion with `FixedWeight(last(...))`.

## Budgets and backends

`NoStop()` disables tolerance acceptance, but still respects the configured budget,
nonfinite detection, and exact completion. With `K≥N`, finite completion requires
`N` fine batches and `N-1` corrections of the stored starts.

Bulk backends support all these criteria. The pipelined backend supports only
fixed-weight `WeightedDefect(p=Inf)`, explicit `FrontierRelativeUpdate`, and
`NoStop()`. It rejects unsupported criteria instead of substituting another test.

Custom criteria subtype `AbstractCriterion` and implement
`NSDETimeParallel.residual(criterion, cache, k)` and
`NSDETimeParallel.accept(criterion, value)`. Also define
`NSDETimeParallel.threshold(criterion)` for reporting, or keep a field named `ϵ`
for the default threshold. See [Extensions](extensions.md) for a tested example.
The core loop does not select concrete
criterion types. A custom adaptive weight rule implements the corresponding
`weight_rule`, `WeightState`, `needs_fine_history`, and `update_weight!` hooks.
