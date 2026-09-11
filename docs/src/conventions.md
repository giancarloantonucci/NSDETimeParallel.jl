# Conventions

Every rule on this page exists because its violation produced a real bug in this ecosystem's history; this is the page to read before touching weights, criteria, or boundaries.

## The weight convention (read this one twice)

There are two natural bases for the exponential weights, and both live in this ecosystem. Confusing them has caused three independent bugs, which is why this box exists.

!!! warning "Per-time vs per-chunk"
    **`Weights.w` is a per-*time* base**: `w = exp(λ)` with λ the amplification rate per unit time. The proximity function applies it through time exponents, `w^(T₁ − Tₙ)`, so the discount a boundary receives is `exp(-λ · elapsed time)` regardless of how time is chunked. The thesis's per-*chunk* base is `w^ΔT = exp(λΔT)`; the two coincide **only when the chunk length is 1**.

    - `Weights(w = exp(λ))`, and the rate measured by `Weights(updatew = true)`, are per-time.
    - [`Wnorm`](@ref) takes the **per-chunk** base (its exponents are chunk indices, matching the thesis's `Wₙ = w⁻ⁿ`); its docstring carries the same warning.
    - Weights are anchored at the **window start** (`T₁`), never at absolute `t = 0`: a window beginning at `t₀ ≠ 0` must not change any discount.

If you are porting a script older than this page: `exp(2^(1/29)·ΔT)` and friends are pre-convention hacks; the ecosystem-wide Lorenz value is `exp(0.9056)` per unit time.

## What the two criteria actually accept

`ψ₁` compares successive iterates; `ψ₂` compares chunk starts against the fine values they should equal, discounted by the weights. They do not merely differ in speed — they accept **different objects**, and at the wrong chunk length each is degenerate in its own way. (For *choosing* between them, and for getting the weight, see the [criterion](criterion.md) page; this section is about what the choice commits you to.)

!!! note "The two-pole degeneracy"
    With chunks much longer than one Lyapunov time of the dynamics:

    - under **`ψ₁`** no meaningful contraction happens until finite termination sweeps through, so runs "converge" at `K ≈ N` — at which point the result *is* the chunked fine solve, computed at roughly `N/2` sweeps of work per solve. You have paid a large multiple of the serial cost to reproduce the serial answer.
    - under **`ψ₂`** the discount `exp(-λ(Tₙ - T₁))` annihilates late boundaries, so runs accept at `K ≈ 1–2` with O(1) discontinuities at late chunk seams: the result is a **pseudo-orbit**, statistically useful for chaotic problems, pointwise wrong.

    The regime in which the criteria mean something is chunks of roughly **one Lyapunov time or less**: `K` small but plural, seams present but bounded.

Two consequences you should design around rather than discover:

- **`ψ₂`'s tail is unconstrained by design.** The discount forgives late boundaries; do not assume they are small — *measure* them with [`seams`](@ref) / [`maxseam`](@ref), and report them next to any statistical claim.
- **Finite termination is exact — in floating point too.** After `k` corrections the first `k` boundary starts are exact, so at `K = N` the result is the chunked fine solve. Two implementation details make that hold bit for bit rather than "up to roundoff": the correction is evaluated as `F + (Gnew − Gold)`, which returns `F` exactly when the coarse values cancel exactly (the other association, `(Gnew + F) − Gold`, rounds `F` into `Gnew` first and loses low bits — an earlier version did that, and the suite carried a `< 1e-13` seam bound with a comment blaming floating point), and the newly exact interface is *assigned* from the fine output instead of recomputed. A fixed-step chunk solve whose step divides the chunk ends on the chunk boundary bit for bit (the last scheduled step is placed on `tN`; see `NSDERungeKutta.fixedstep_count`), so the boundary lookup returns the final node rather than interpolating across a one-ulp overshoot. The suite asserts `maxseam == 0.0` at `ϵ = 0, K = N`. Backends agree bitwise with each other as before. [`contractionrate`](@ref) still excludes zeros and non-finite entries: an exact zero at the fixed point is a property of termination, not a rate.

## Boundary ownership: `flatten`

Adjacent chunks share a boundary instant with two values at it: chunk `n`'s **fine endpoint** and chunk `n + 1`'s **corrected start** (they differ by exactly the seam). [`flatten`](@ref) keeps the **fine-endpoint side** — chunk 1 in full, then `u[2:end]` of every later chunk, so the dropped value is each later chunk's corrected start — producing one point per time on the fine grid. If you concatenate chunks by hand you will either duplicate boundaries or make the opposite ownership choice silently; the historical `collapse` helpers did both *in the same script*. Use `flatten`, and audit what it hid with [`seams`](@ref).

## Measure, don't guess

Every quantity the cost and tracking formulas consume has a measuring function; stage-count arithmetic and step-ratio guesses are how the old figures went wrong.

| Quantity | Meaning | Measure with |
|---|---|---|
| ζ | coarse/fine cost ratio per chunk | [`costratio`](@ref) |
| β | per-iteration error contraction | [`contractionrate`](@ref) |
| S | predicted speedup at K iterations | [`theoretical_speedup`](@ref) |
| K̂ | largest K meeting a speedup target | [`iteration_budget`](@ref) |
| seams | what was actually accepted | [`seams`](@ref), [`maxseam`](@ref) |

Match the ceiling to the measurement: [`theoretical_speedup`](@ref)'s default is the *work-credited* bound (the shrinking `k:N` sweeps are saved **work**); a wall-clock benchmark with a full worker pool does **not** collect that credit — a parallel sweep costs one chunk-time however few chunks remain — so judge timings against `estimate = :wallclock`, and judge the pipelined backend against `estimate = :aubanel`. Judging a measurement against the wrong ceiling flatters or slanders it; the [backends](backends.md) page spells out which is which.

## Execution modes

The primitive-based backends change wall time only: results are bitwise identical across `SERIAL`, `THREADS`, `DISTRIBUTED` and `MPI`, and the test suite enforces it. `PipelinedMPIBackend` is the documented exception — the same fixed point reached under frontier convergence semantics. Operational guidance (which backend, when, and the traps) lives on the [backends](backends.md) page.

## Plotting

Plotting goes through Plots (RecipesBase); attributes are `variables`, `iscomplex` and `skip`. Parallel and windowed solutions draw chunk-by-chunk and window-by-window — seams and overlaps are shown, not hidden, per this page's second section.
