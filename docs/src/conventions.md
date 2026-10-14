# Representation and accounting

## Starts, outputs, and reference

Parareal stores `N` chunk starts. Each returned chunk is a fine solve from its
stored start; the final time value comes from the last chunk. `boundaryvalues`
returns the starts plus that terminal output. `seams` measures jumps between a
chunk's fine output and the next stored start.

The chunk span is `(tN-t0)/N`. NSDERungeKutta prepares one nominal schedule:
full steps plus a shortened last step when necessary. All translated chunks use
that schedule. Moving windows pass their nominal duration explicitly when resizing;
translation alone does not rebuild caches. Chunk output is read from `u[end]`; generic non-NSDE solver adapters
retain interpolation only if their returned trajectory overshoots the boundary.

`reference(problem, parareal)` applies the fine solver serially on this same grid.
For deterministic finite solves, exact completion reproduces these interface values
bitwise. The correction forms `F + (Gnew-Gold)` and directly assigns the newly exact
start. The guarantee concerns the discrete reference, not the exact ODE trajectory.

## Histories and reusable configuration

`history=:none` is the default. `:interfaces` saves boundary vectors for each
batch in `diagnostics.interfaces`; `:iterates` saves full chunk trajectories in
`solution.iterates`. Residual and weight histories are always in diagnostics.
`errors` aliases the residuals in diagnostics. Use `history` to choose stored output.

`Parareal`, its parameters, criterion, and weight rules are immutable configurations.
Cached states and solver workspaces belong to a run. `remake(parareal; criterion=...,
finesolver=..., coarsesolver=..., N=..., K=...)` changes selected settings.
A cache can switch weight-rule and fine/coarse solver types, including fixed and
adaptive RK solvers; method-dependent storage is rebuilt when needed. The chunk
count and state type/shape must stay fixed. Rebuild the cache when those change,
and create a new solution buffer if the new solver uses a different output type. Running a solver does not mutate its settings.

## Counts and status

For bulk execution, `K` limits fine batches. Every finite return has
`fine_batches=sweeps+1`, for both relative-update and weighted-defect criteria.
Exact completion uses `N` batches and `N-1` corrections. The final batch evaluates
the starts that are actually returned; no untested final correction is made.

Diagnostics record actual fine/coarse chunk calls after exact-prefix reuse.
Statuses are `:converged`, `:exact`, `:budget_exhausted`, and `:nonfinite`.
A numerical `SolverFailure` (including `NewtonFailure` or `StepFailure` from
exhausted adaptive retries) is logged with its time span and backtrace, and the
failed chunk is replaced by NaNs. All backends finish with `:nonfinite`.
Other exceptions, including programming and configuration errors, propagate on
serial, threaded and distributed backends (possibly wrapped by the backend).
Under MPI, all ordinary chunk exceptions become logged `:nonfinite` returns so
peers can finish their communication instead of hanging. This covers the coarse
initial pass and corrections as well as fine solves. Standalone Runge–Kutta
solves still throw numerical failures. Interrupts, stack overflow and out-of-memory
failures are rethrown; an active multi-rank MPI backend aborts for these fatal
failures. Nonfinite results can
terminate partway through a computation, so the finite-return count identity
should not be used to reconstruct their work. If the coarse pass fails before
any fine batch, the residual and measured-weight histories are empty. Unsolved
chunk trajectories carry NaNs on their proper time spans, so interpolation and
`boundaryvalues` cannot expose uninitialised or stale output. If a correction
fails, only the already exact prefix is retained.

Pipelined runs have no global batches: `sweeps` and `fine_batches` are `nothing`.
Use `local_sweeps` and the actual total call counts instead.

## Work models

`theoretical_speedup(K, N, ζ)` takes **fine batches** K. A finite bulk return has
s=K−1 corrections. With the exact-prefix optimization, the coarse call count is
`K_G=K*(N-1)-K*(K-1)/2`. The default `estimate=:wallclock` gives
`N/(K+ζ*K_G)`, charging one fine-chunk time per batch and all sequential coarse
calls. It assumes enough workers and equal chunk costs, and excludes communication,
synchronization, and load imbalance; it is not an unconditional timing ceiling.

For a measured finite run, prefer:

```julia
ζ = costratio(problem, parareal)
S = theoretical_speedup(solution.diagnostics, ζ)
# Equivalently: theoretical_speedup(solution, ζ)
K = iteration_budget(2.0, numchunks(solution), ζ)
```

The diagnostics form uses recorded calls directly, including those from a
zero-correction acceptance or exact completion. `iteration_budget` returns the
largest fine-batch budget meeting the requested modelled speed-up; its correction
budget is K−1, and K=0 means that the speed-up target is infeasible. The helper
says nothing about whether that budget is sufficient for convergence.

`estimate=:work` divides the actual fine-call count by N instead of charging a
whole batch. It assumes perfect redistribution of the shrinking workload, which
the bulk backend does not provide. Use the same model for both helpers. The old
budget names `:naive` and `:improved` alias `:wallclock` and `:work`, now with the
corrected coarse counts. The speed-up helper also retains `wallclock=false` as
an explicit spelling of `:work`.

Nonfinite, unstarted, and pipelined diagnostics are refused. Pipeline calls do not
have a global batch count, and total calls alone cannot give their critical path.
The old `:aubanel` expression has been removed from this helper rather than
presented as a bound for the native pipeline.
