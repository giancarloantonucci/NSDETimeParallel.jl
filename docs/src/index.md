# NSDETimeParallel.jl

Parareal with explicit stopping targets, chunk schedules, and work diagnostics.

## First example

```julia
using NSDERungeKutta, NSDETimeParallel
problem = Logistic(0.3, (0.0, 2.0))
par = Parareal(RK4(h=0.005), RK4(h=0.05); N=4)
result = solve(problem, par; backend=ThreadsBackend())
@assert result.diagnostics.status in (:converged, :exact)
result.diagnostics.fine_batches
```

Start with the [guided tour](guide.md), then consult:

- [Criteria](criterion.md)
- [Conventions and cost models](conventions.md)
- [Backends](backends.md)
- [Extensions](extensions.md)
- [Quick reference](quick_reference.md)
- [API](api.md)

## Environment

These pages describe the current development checkout (Julia 1.10 or later).
Use compatible versions of the NSDE packages in one project. For local
checkouts, develop each package with `Pkg.develop(path="/path/to/package")`;
`Pkg.status()` reports what will actually load. A registered release may expose
an older API. Plotting examples additionally need Plots.

## Reference
G. A. Antonucci (2025), *Solving chaotic ODEs with time-parallel algorithms*,
DPhil thesis, University of Oxford. DOI: [10.5287/ora-xoevp0dez](https://doi.org/10.5287/ora-xoevp0dez).
The thesis provides the mathematical background; the API and implementation
conventions are specified by these package docs.
