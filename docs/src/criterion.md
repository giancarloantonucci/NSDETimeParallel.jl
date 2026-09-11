# Choosing the convergence criterion

Parareal stops when `ψ(iterate) ≤ ϵ`. Both live in [`Tolerance`](@ref), and the choice of `ψ` is the most consequential user-facing decision in this package: it decides what "converged" *means*, and therefore how many iterations you pay for and what errors you silently accept. This page is the practical guide; the docstrings hold the formulas.

## ψ₁ — the honest default

[`ψ₁`](@ref) is the mean relative boundary update, `mean(‖U[n] − V[n]‖ / ‖U[n]‖)` against the previous iterate. It is uniform — every chunk boundary must settle — and needs no tuning. Use it whenever the time span is short relative to the dynamics' growth rate, and always when you need "converged" to mean *converged everywhere*.

Its failure mode is long chaotic spans. When perturbations grow like `e^{Λt}`, a boundary far into the window cannot settle to ϵ while its upstream boundaries are still moving by more than `ϵ·e^{−Λ(Tₙ − T₁)}` — the criterion demands the impossible, and the iteration count is driven toward the finite-termination diagonal `K = N`, at which point Parareal is a very elaborate serial solve.

## ψ₂ — the Lyapunov-discounted criterion

[`ψ₂`](@ref) weights boundary `n`'s error by `w^{T[1] − T[n]}`. With `w = exp(Λ)` for a problem with Lyapunov rate Λ, the discount matches the amplification: an error the dynamics will grow by `e^{Λ(Tₙ − T₁)}` is discounted by exactly that factor, so every boundary is judged by what its error *means* propagated to a common reference, not by its raw size. This is the moving-window criterion — it is what lets NSDEMovingWindow accept a window whose early stretch (the part the next window will keep) has settled.

State the trade honestly: convergence in ψ₂ means *early boundaries settled*, not *uniformly within ϵ*. Late-window seams are unconstrained **by design**. Whenever a ψ₂-converged result feeds a claim about the whole span, audit what was forgiven: [`seams`](@ref) and [`maxseam`](@ref) report the boundary discontinuities directly, and the true terminal error against a serial fine solve is the final judge.

ψ₂ also has a variational reading worth knowing: it is (up to squaring) a weighted least-squares objective on the boundary residuals, `φ(U) = ½Σ‖e^{−λ(Tₙ−T₁)}(Uₙ − F(Uₙ₋₁))‖²` — Parareal acts as an iteration driving that objective down, which is where the criterion historically comes from (an optimisation-formulation study that also compared Gauss–Newton and quasi-Newton methods on the same objective).

Two mechanical notes: `w` is a **per-unit-time base** (`w = exp(Λ)` for per-time rate Λ — the exponent is a time difference, so the base must be per-time; the [conventions](conventions.md) page has the per-time/per-chunk warning box), and ψ₂ clamps `w` to at least 1: a sub-unit base would *inflate* late errors, which is never the intent.

## ψ∞ — the local (pipelined) criterion

[`ψ∞`](@ref) is ψ₂ with the mean replaced by the **maximum**: the largest weighted boundary defect, same defect, same per-unit-time base, so `ψ₂ ≤ ψ∞ ≤ N·ψ₂` always. Two reasons to reach for it (thesis §3.4, `subsec:local_proximity`):

- **Sharper certification.** `ψ∞ ≤ ϵ` certifies the weighted-max error within `s(θ_F)·ϵ`, where `θ_F = Λ_F/w` and `s(θ_F) = Σ_{j<N} θ_F^j` (bounded by `1/(1 − θ_F)` when `w > Λ_F`) — the ℓ¹ constant carries an extra factor `N`, so per-chunk certification does not loosen with the processor count. To certify a target radius `r`, set `ϵ = r/s(θ_F)`.
- **It is the one criterion a pipeline can honour.** Acceptance under a max decomposes into per-chunk tests conjoined left to right — one boolean riding the messages `PipelinedMPIBackend` already sends. A mean is a shared budget: no chunk can accept alone, so enforcing ψ₂ would re-serialise the pipeline.

The trade is the mirror of ψ₂'s: a max cannot average away one stubborn boundary, so on bulk-synchronous backends ψ∞ is the *stricter* stopping rule at equal ϵ. Late-window seams remain discounted by design, exactly as under ψ₂ — audit with [`seams`](@ref)/[`maxseam`](@ref) as ever.

## Getting the weight

In order of preference:

1. **A known rate.** For textbook systems Λ is literature (Lorenz: Λ ≈ 0.9056, so `Weights(w = exp(0.9056))`).
2. **A probe.** Solve one representative window with `Weights(updatew = true)`, read off the measured base, and **freeze** it:

```julia
probe = Parareal(finesolver, coarsesolver; parameters = params,
                 tolerance = Tolerance(ϵ = ϵ, ψ = ψ₂, weights = Weights(updatew = true)))
solve(copy(problem, u0, t0, t0 + τ), probe)
ŵ = probe.tolerance.weights.w                    # measured on iterate data
frozen = Tolerance(ϵ = ϵ, ψ = ψ₂, weights = Weights(w = ŵ))
```

3. **Not** `updatew = true` left on in production. The updater floors `w` at the largest rate it has measured, so it can only ratchet the discount up — and under NSDEMovingWindow's Zoom strategy it silently undoes every tightening, which is why the `MoWi` constructor refuses that pairing outright.

A characterisation worth knowing (pinned in the test suite): on *smooth orbit data* the updater's quotient spans two chunk gains and returns the **square** of the one-chunk gain. The probe regime is different — mid-iteration Parareal data is iterate discrepancies, not orbit spacings — but treat a probe-measured `ŵ` as an estimate to sanity-check against what you know of the dynamics, not as gospel.

## What the criterion cannot see

The K you pay is set by the criterion; whether that K is *worth paying* is set by the speed-up ceiling. [`iteration_budget`](@ref)`(S, N, ζ)` inverts the efficiency estimates: the largest K that still meets a target speed-up S. If your criterion routinely converges above that budget, no backend will save the run — loosen ϵ, shorten the window, or accept a smaller S.

Finally, backend interaction: `PipelinedMPIBackend` honours exactly one ψ — [`ψ∞`](@ref), whose prefix-decomposable acceptance rides the pipeline's own messages (with a **preset** weight; the updater needs a global view and is rejected there). Any other ψ is ignored in favour of the per-boundary stagnation frontier (see the Backends page), because a shared-budget criterion evaluated every sweep would re-serialise the very pipeline it runs on.
