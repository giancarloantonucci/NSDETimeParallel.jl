# Extending Parareal

Start with a working serial solver and the built-in criteria. Extensions use
qualified method names so that the package whose interface is extended is clear.
The executable companion is `examples/extensions.jl`.

## A custom maximum-seam criterion

```@example extensions
using NSDERungeKutta, NSDETimeParallel, LinearAlgebra
struct MaxJump <: AbstractCriterion
    tolerance::Float64
end
function NSDETimeParallel.residual(c::MaxJump, cache, k)
    return maximum(norm(cache.U[n] - cache.F[n]) for n in eachindex(cache.U))
end
NSDETimeParallel.threshold(c::MaxJump) = c.tolerance
NSDETimeParallel.accept(c::MaxJump, value) = isfinite(value) && value <= c.tolerance
problem = Logistic(0.3, (0.0, 2.0))
p = Parareal(RK4(h=0.005), RK4(h=0.05); N=4, criterion=MaxJump(1e-8))
r = solve(problem, p)
@assert r.diagnostics.status in (:converged, :exact)
nothing # hide
```

`U[n]` is start `n`; `F[n]` is the fine endpoint from chunk `n-1`, with
`F[1]=U[1]`. Thus the first defect is zero. The loop also records the criterion threshold. Define `threshold` and
`accept` explicitly, as here, or supply a field named `ϵ` and use their defaults. This criterion has the
same unweighted maximum-defect target as `WeightedDefect(p=Inf,weight=FixedWeight(1))`;
it illustrates the interface rather than adding a new numerical method.

Custom criteria work on bulk backends. A pipeline criterion needs a separate
protocol; implementing `residual` does not make it pipeline-compatible. Cache
fields used by an extension are implementation-facing and should be tested when
upgrading the package.

## Wrap an existing serial integrator

```@example extensions
struct Wrapped{S} <: AbstractInitialValueSolver
    inner::S
end
NSDEBase.initialize_cache(problem::AbstractInitialValueProblem, s::Wrapped) =
    initialize_cache(problem, s.inner)
NSDEBase.initialize_solution(problem::AbstractInitialValueProblem, s::Wrapped; kw...) =
    initialize_solution(problem, s.inner; kw...)
NSDEBase.solve!(cache, sol, problem::AbstractInitialValueProblem, s::Wrapped) =
    solve!(cache, sol, problem, s.inner)
NSDEBase.solve(problem::AbstractInitialValueProblem, s::Wrapped; kw...) =
    solve(problem, s.inner; kw...)
NSDEBase.chunk_solver(s::Wrapped, duration::Real) = Wrapped(NSDEBase.chunk_solver(s.inner, duration))
pw = Parareal(Wrapped(RK4(h=0.007)), RK4(h=0.05); N=4, criterion=NoStop())
rw = solve(problem, pw)
@assert boundaryvalues(rw.lastiterate) == reference(problem, pw)
nothing # hide
```

This wrapper delegates scheduling as well as stepping. A foreign integrator
needs an adapter translating its cache, solution, cached `solve!`, and one-off `solve` operations, not just
a new type name. Its solution must expose writable time/state storage and the
indexing/interpolation operations the Parareal solution expects. Check early
acceptance, non-dividing chunk lengths, exact completion, and failure behaviour.
Scheduled endpoints should be reached directly; overshooting adapters may use
interpolation and should not assume bitwise reference agreement automatically.

Throw a `SolverFailure` subtype for numerical failure. Leave programming errors
as ordinary exceptions. Bulk non-MPI backends propagate ordinary errors; MPI
logs ordinary chunk exceptions and returns `:nonfinite` to let peers complete
communication. See [Conventions](conventions.md) for fatal errors and output invalidation.

For a custom problem, provide writable `rhs`, `u0`, `tspan` and an alias-preserving
`remake(problem; u0, tspan)` as described in NSDEBase's guide. A custom solver for
Stretch additionally needs the step-resizing interface when adaptation uses it.
