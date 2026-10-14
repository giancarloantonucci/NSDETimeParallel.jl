# Glossary and quick reference

| Term | Meaning in a bulk Parareal run |
|---|---|
| Chunk | One of `N` time intervals |
| Start | A stored initial value for a chunk; there are `N` |
| Boundary values | The starts plus the terminal fine output; there are `N+1` |
| Fine batch | A round of fine propagation from the current starts |
| Sweep | A coarse correction between fine batches |
| Seam | A fine endpoint minus the neighbouring stored start, up to sign |
| Exact completion | Matching the chunked fine reference after at most `N` batches |
| Reference | `reference(problem, parareal)`, not the exact ODE solution |
| Weight | Per-unit-time discount base; diagnostics store the value used by each batch |
| `ζ` | Measured coarse/fine cost ratio for one chunk |
| Window | A MoWi subproblem, possibly overlapping its predecessor |
| Clipped | A genuinely shortened endpoint window, excluded from default MoWi statistics |

For finite bulk returns, `fine_batches=sweeps+1`. Partial failures need actual
recorded counts. A pipeline has local iterations and no global batch count.

## Common calls

The following is a reference list; identifiers such as `problem` and `par` must
first be constructed as in the [guided tour](guide.md).

| Task | Call |
|---|---|
| Default mixed update test | `RelativeUpdate(atol=1e-12, rtol=1e-9)` |
| Mean discounted defect | `WeightedDefect(ϵ=1e-9, weight=FixedWeight(1))` |
| Maximum discounted defect | `WeightedDefect(ϵ=1e-9, weight=FixedWeight(1), p=Inf)` |
| Empirical same-chunk weight | `SecantWeight(w=1, τ=1e-10)` |
| Disable early acceptance | `NoStop()` (budget still applies) |
| Save boundary histories | `solve(problem, par; history=:interfaces)` |
| Save full iterates | `solve(problem, par; history=:iterates)` |
| Inspect continuity | `seams(r)`, `maxseam(r)` |
| Compare to reference | `Wnorm(r, reference(problem, par))` |
| Read measured weight | `last(r.diagnostics.weights)` after checking it is nonempty |
| Estimate cost | `theoretical_speedup(r, costratio(problem, par))` |
| Change configuration | `remake(par; N=8, criterion=NoStop())` |

`Wnorm` measures a weighted interface distance using the last measured base.
Its reference solve is extra work for validation, not part of ordinary stopping.

## Diagnose before tuning

| Observation | Check |
|---|---|
| Initial `Inf` relative-update residual | Expected before a corrected previous iterate exists; read status |
| Empty residual history | The coarse pass may have failed before a fine batch |
| Batch count near `N` | Coarse stability/accuracy, chunk length, phase errors, and tolerance |
| Small count but large seams | Averaging or discounting may permit this; inspect unweighted defects |
| Model speed-up below one | Coarse cost can outweigh saved batches |
| `contractionrate` is `NaN` | Fewer than two positive finite residuals are available |
| Full-history request fails on reuse | The solution was not allocated with `history=:iterates` |
| Cache reuse fails | Check state type/shape, chunk count, time span, and solver output type |

## Reproducible and meaningful measurements

Record package versions and source revisions, Julia version, platform, backend,
thread/rank counts, problem parameters, criteria, and floating-point types.
Retain status, residuals, weights, fine/coarse call counts, and accuracy diagnostics.
Save seeds whenever initial conditions, forcing, or analysis use randomness.
Deterministic tests can establish agreement for a specified configuration; they
do not promise bitwise agreement across platforms, BLAS libraries, compiler
versions, random right-hand sides, or different operation orders.

Store portable tables alongside Julia-specific serialized records. Only deserialize
trusted Julia files. Include a reference-completion check and a step-refinement
check appropriate to the accuracy claim. Benchmark warmed-up runs and state
whether setup, collection, and communication are included. Cache reuse reduces
setup cost but does not make it zero. The tiny tutorial problems teach the API;
they are not demonstrations of parallel speed-up.
