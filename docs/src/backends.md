# Execution backends

Select with `backend=...`. The string-valued `mode` keyword has been removed.

| Backend | Scheduling | Criteria | History |
|---|---|---|---|
| `SerialBackend()` | bulk, calling process | all bulk rules | all |
| `ThreadsBackend()` | bulk, threads | all bulk rules | all |
| `DistributedBackend()` | bulk, worker pool | all bulk rules | all |
| `MPIBackend()` | bulk, MPI ranks | all bulk rules | all |
| `PipelinedMPIBackend()` | asynchronous frontier | fixed maximum defect, frontier update, no stop | residual/local-count diagnostics |

Start Julia with `--threads=N` for threads. Adaptive NSDERungeKutta solves keep
the current step size and history local to each solution, so they do not mutate a
shared chunk solver. User right-hand sides must likewise avoid shared mutable state.

For Distributed, start workers with the same project and load the packages on each:

```julia
using Distributed
addprocs(2; exeflags="--project=.")
@everywhere using NSDERungeKutta, NSDETimeParallel
solution = solve(problem, parareal; backend=DistributedBackend())
```

MPI requires `MPI.Init()` and an SPMD launch. MPI remains a required dependency
in this release. See `test/mpi/runtests_mpi.jl` for an executable example.

## Pipeline contract

The pipeline overlaps coarse propagation and fine work, passing a convergence
frontier from left to right. It does not have a shared current iterate at every
instant. Consequently `RelativeUpdate`, `MeanRelativeUpdate`, averaged defects,
and changing secant weights are rejected rather than silently replaced.

Use `WeightedDefect(ϵ=..., p=Inf, weight=FixedWeight(...))` for maximum-defect
acceptance. A rank waits for the incoming start to freeze before accepting its
outgoing defect. `FrontierRelativeUpdate(ϵ=...)` explicitly selects the legacy
relative outgoing-update heuristic. `NoStop()` runs to the budget or diagonal.
At least `N` ranks are required; surplus ranks take part only in collectives.

`diagnostics.local_sweeps` records each chunk's local iteration count. Global
fine/coarse call totals include final fine solves needed after a late start update.
`fine_batches` and `sweeps` are `nothing`. Residuals aggregate ranks active at each
local iteration index and must not be interpreted as bulk residuals.

Full or interface histories would require additional synchronization and are
rejected. With `NoStop()` and `K≥N`, the final collected result matches the chunked
fine reference. For performance comparisons use the matching scheduling model.

Under `MPIBackend`, all ranks receive a residual history trimmed to the completed
fine-batch count, and `solution.errors === solution.diagnostics.residuals` holds
on every rank, including zero-batch failures and reused solutions.

The pipelined backend initializes interior starts through its coarse communication
chain. It does not use `makeGs` overlap guesses from a window driver; only the
window's initial value is preserved. Use a bulk backend for interior warm starts.


## Complete launch examples

`examples/distributed.jl` starts two workers in the active project, loads NSDE on
each, checks the chunked fine reference, and removes the workers. Load named
right-hand-side functions on every worker too. Each worker needs access to the
same package environment and any data used by its right-hand side.

`examples/mpi.jl` runs the bulk backend with one rank per chunk and checks the
returned histories on every rank. With MPI.jl's compatible launcher installed:

```sh
mpiexecjl -n 2 julia --project examples/mpi.jl
```

Launch all ranks into the same project. Both MPI backends require at least `N`
ranks. Do not run a multi-rank program by starting a single Julia process and
expecting `MPIBackend()` to create workers. The pipeline is a separate protocol,
not an interchangeable performance switch; see its supported criteria above.

For threaded matrix-heavy workloads, consider `BLAS.set_num_threads(1)` after
`using LinearAlgebra`, then measure: nested BLAS and chunk parallelism can
oversubscribe cores. Each chunk must own any mutable scratch storage used by its
right-hand side. Thread safety is a property of user functions as well as the solver.
