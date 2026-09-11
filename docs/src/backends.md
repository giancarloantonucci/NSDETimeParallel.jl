# Execution backends

The Parareal algorithm is written once (`src/parareal/parareal.jl`); a backend supplies only execution primitives, so a fix lands in all four primitive-based backends at once and they cannot drift. Pick one with `mode = "SERIAL" | "THREADS" | "DISTRIBUTED" | "MPI" | "PIPELINED"` or by passing a backend object. The fifth backend, `PipelinedMPIBackend`, is deliberately **not** an adapter but a documented semantic variant with its own loop — see its section below.

## SerialBackend (default)

Everything in the calling task. The reference implementation: the test suite asserts the other backends agree with it **bit for bit**. Also the right choice when the fine solve is too cheap to parallelize.

## ThreadsBackend

Fine chunks over `Threads.@threads :dynamic` in shared memory — the first thing to try on one machine. Start Julia with `--threads=N`.

## DistributedBackend

Fine chunks over `pmap` on the default worker pool (Distributed.jl). The schedule adapts to however many workers exist, so the result does not depend on the worker count. Chunk solutions come back through memory, not disk. Load the packages on the workers first:

```julia
using Distributed
addprocs(4)
@everywhere using NSDETimeParallel, NSDERungeKutta
solution = solve(problem, parareal; mode = "DISTRIBUTED")
```

## MPIBackend

SPMD over MPI ranks, one rank per chunk: boundary values through point-to-point messages, final chunks through a generic gather — all in memory. Needs `MPI.Init()` and `mpiexec`; `test/mpi/runtests_mpi.jl` is a ready-made recipe:

```
mpiexecjl -n 4 julia --project=@. test/mpi/runtests_mpi.jl
```

Disk is used only for the opt-in `saveiterates` history under MPI, in a per-run temporary directory by default.

## PipelinedMPIBackend — a semantic variant

Pipelined (task-scheduled) Parareal over MPI, one rank per chunk, after Aubanel's scheduling (thesis §2.7): there is **no root**. The serial coarse chain is passed rank to rank, and every rank runs its next fine solve *before* blocking on the corrected start from its left neighbour — the fine work hides the serial coarse cost, which is exactly what the root-centric `MPIBackend` cannot do.

The price is stated openly rather than papered over:

- **Frontier convergence, two flavours.** A chunk stops once every chunk to its left has accepted *and* its own per-chunk test passes, or it reaches the finite-termination diagonal. With `ψ = ψ∞` the test is the **certified** weighted defect of thesis §3.4 (`w^(T[1]−T[n+1])‖uold − F(start)‖ ≤ ϵ`, evaluated with the thesis's one-stage lag; acceptance is the prefix conjunction `A_n = A_{n−1} ∧ t_n`, riding the existing messages) — a converged run satisfies `ψ∞ ≤ ϵ`, hence a weighted-max error within `s(θ_F)·ϵ`. The weight must be **preset** (`updatew = true` is rejected: no rank sees all boundaries at a common sweep). With any other ψ the test is relative stagnation of the outgoing boundary, the global ψ is **ignored** (warning unless `ψ₁`) — a shared-budget criterion needs every boundary once per sweep, which re-serialises the pipeline.
- `solution.errors[k]` records the largest per-chunk test value among chunks still active at sweep `k` — the running restriction of ψ∞ to active chunks, or the largest relative update.
- `saveiterates` is unsupported: the history would force per-sweep synchronisation.

What is preserved: accepted chunks are fine solves from their accepted starts, so `flatten`, `seams` and all downstream statistics behave identically; and at `ϵ = 0`, `K = N` finite termination makes the result the chunked fine solve, bitwise — the same fixed point as `SerialBackend`, reached by a different route (`test/mpi/runtests_mpi.jl` asserts this).

```
mpiexecjl -n 4 julia --project=@. test/mpi/runtests_mpi.jl
```

Use it when the coarse solver is expensive enough that the root's serial chain dominates (`ζN` large); judge its timings against `theoretical_speedup(K, N, ζ; estimate = :aubanel)`, never against the root-centric ceilings.

## Choosing a backend in practice

Everything below is measured, not folklore.

- **The guarantee first: the primitive-based backends change wall time only.** Results are bitwise identical across `SERIAL`, `THREADS`, `DISTRIBUTED` and `MPI` — the test suite enforces it. (`PipelinedMPIBackend` is the documented exception: same fixed point, different convergence semantics — see its section.) Use `SERIAL` for anything whose point is determinism or diagnostics; you lose nothing but time.
- **Chunk work must dwarf scheduling.** Every sweep pays a fixed per-dispatch cost; if a chunk's fine solve is microseconds, the `THREADS` task spawn — let alone `DISTRIBUTED`'s serialisation — exceeds the work and parallelism slows you down. Measure your machine's dispatch floor with a zero-work probe (a Parareal solve whose fine solver takes one trivial step) before blaming the algorithm.
- **What limits Parareal is not the barrier.** An embarrassingly parallel workload under the same fork–join machinery scales near-ideally; Parareal's cap is its serial coarse fraction — the `ζN/(1 + ζN)` of the cost model — plus dispatch overhead on small chunks.
- **`MPIBackend` is SPMD** — one rank per chunk, `MPI.Init()` and `mpiexec` required. It is the wrong tool for interactive work and plotting sweeps; for worker pools over MPI transport driven from a script, use MPIClusterManagers' `MPIWorkerManager` with the `DISTRIBUTED` backend instead — and mind the id mapping: Distributed worker 2 is MPI rank 0.
- **Pool set-up details that bite:** `addprocs(...; topology = :master_worker)` skips the all-to-all connection set-up the `DISTRIBUTED` backend never uses (noticeable at large pools). On heterogeneous CPUs (Apple silicon P/E cores) a "full" thread pool mixes fast and slow workers, and the slowest chunk sets the sweep time. Physical cores (`Hwloc.num_physical_cores()`), not `Sys.CPU_THREADS`, are the honest pool size for compute-bound chunks.

## Judging speed-ups

Parareal has its own ceiling: [`theoretical_speedup`](@ref)`(K, N, ζ)` for `K` iterations over `N` chunks with a **measured** coarse/fine cost ratio ζ — [`costratio`](@ref) measures it for you on one chunk. Three estimates are available: the work bound (default), the wall-clock bound (`estimate = :wallclock`, for timing benchmarks with a full worker pool), and the Aubanel bound (`estimate = :aubanel`, the comparator for `PipelinedMPIBackend`). Judge real timings against the matching ceiling, never against the serial fine solve wishfully. `benchmark/scaling.jl` and `benchmark/scaling_sweep.jl` produce pool-scaling measurements (speed-up, and efficiency across thread counts); they print measurements, not bounds — the ceilings above are what those curves should be judged against.
