# NSDETimeParallel.jl

Parareal with explicit stopping targets, chunk schedules, and work diagnostics.



## Start here

```julia
using NSDERungeKutta, NSDETimeParallel
problem = Logistic(0.3, (0.0, 2.0))
par = Parareal(RK4(h=0.005), RK4(h=0.05); N=4)
result = solve(problem, par; backend=ThreadsBackend())
@assert result.diagnostics.status in (:converged, :exact)
result.diagnostics.fine_batches
```

Start Julia with `--threads=4` for this example; `SerialBackend()` is the default.
Parareal constructs the fine/coarse chunk schedules itself. The default
`RelativeUpdate()` requires every stored start to pass a mixed absolute/relative
update test. `WeightedDefect` deliberately changes the acceptance target.

For deterministic finite bulk solves, `NoStop()` with `K=N` reaches the matching
chunked fine reference in `N` fine batches and `N-1` corrections. This is discrete
exactness, not exactness relative to the ODE. Use `reference(problem, par)` to
check it. Read status before interpreting residuals or speed-up estimates.

- [Criteria and weights](docs/src/criterion.md): per-unit-time weights and tail limits.
- [Conventions](docs/src/conventions.md): histories, counts, reuse, failures, and cost models.
- [Backends](docs/src/backends.md): serial, threads, Distributed, and MPI.
- [Extensions](docs/src/extensions.md): a complete criterion and solver wrapper.
- [Quick reference](docs/src/quick_reference.md): glossary and common calls.

MPI remains a required dependency. Pipelined MPI has a distinct frontier protocol:
it has no global batch count and supports only a subset of criteria.

## Learn by running

The [guided tour](docs/src/guide.md) explains each step. Its companion
[`examples/guided_tour.jl`](examples/guided_tour.jl) runs the numerical examples
with assertions. Run it in a project containing this package and its example
dependencies; use Julia 1.10 or later.

The ecosystem is NSDEBase → NSDERungeKutta → NSDETimeParallel → NSDEMovingWindow.
ArrowMacros provides supporting macros; they are not required in user examples.

## Installation and versions

For a registered release, use `Pkg.add("NSDETimeParallel")`. These docs describe the
current development API; a registered version may differ. To use a local
checkout, activate your working environment and `Pkg.develop(path="/path/to/NSDETimeParallel")`.
Develop matching local dependency checkouts in that environment too, then check
`Pkg.status()`. The serial and parallel guides also require NSDERungeKutta;
plotting requires Plots.

See the [documentation](https://giancarloantonucci.github.io/NSDETimeParallel.jl/dev/)
and [API reference](docs/src/api.md) for details.

## Build the manual locally

From the package root, with matching dependency checkouts already available:

```sh
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=docs --threads=4 docs/make.jl
```

Documenter executes the guided-tour examples during the build. For local
unreleased dependencies, develop their paths into `docs` before instantiating.
The deployment step only publishes in an appropriately configured CI environment.
