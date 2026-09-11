# Changelog

## 0.2.0

Requires NSDEBase 0.3.1 and NSDERungeKutta 0.2 (test dependency).

### Fixed
- `Wnorm` used the exponent `Ts[n]/(Ts[n−1] − Ts[n])`; it now computes the
  thesis weight `w^{-j}` with `j` window-relative.
- `ψ₁` divided by a zero norm and snapshotted the previous iterate with `.=`,
  aliasing the inner arrays.
- Parareal correction evaluates `F + (Gnew − Gold)` and assigns the interface
  freed by each sweep from the fine output: finite termination is bitwise
  exact and one coarse solve per sweep is saved.
- `update!` skips unresolved ratios and never stores a non-finite weight; the
  per-chunk ratio is converted to a per-time base.
- `iteration_budget` inverts the monotone work bound by bisection on `1:N`,
  clamped to `0:N`.
- `PararealParameters` validates `N, K ≥ 1` and warns on `K > N`.

### Changed (breaking)
- Backends: `SerialBackend`, `ThreadsBackend`, `DistributedBackend`,
  `MPIBackend`, `PipelinedMPIBackend` share one `parareal!`; the `mode=`
  strings still work. The per-mode source files are gone.
- `PararealCache` preallocates chunk problems, per-chunk fine caches and the
  coarse cache/solution; `shiftwindow!` is the only sanctioned way to move a
  cache to a new window.
- `collect!` → `collect_iterates!`, which also copies below-diagonal chunks
  forward; it requires a solution built with `saveiterates = true`.
- `Weights.w` is a finite, positive scalar (vectors had no meaning downstream);
  `δ` is finite and positive and divides the measured candidate only;
  `update!` takes the chunk grid `T`.
- Supported Julia: `1.6` and later (was `1.10`). Verified on 1.6–1.13. The
  threaded fine map uses the default `@threads` schedule: dynamic from 1.8,
  static on 1.6–1.7.

### Added
- `ψ∞`; `boundarytimes`, `boundaryvalues`, `flatten`, `seams`, `maxseam`;
  `theoretical_speedup`, `costratio`, `iteration_budget`, `contractionrate`;
  `NSDEBase.initialize_cache`/`initialize_solution` for `Parareal`.
- Thread-scaling benchmark; MPI test runner (run by a dedicated CI job); docs
  pages on backends, the stopping criterion and conventions; Aqua in the
  test suite.

### Migration
- `collect!(solution)` → `collect_iterates!(solution; directory)`.
- `Weights(w=[...])` → a single scalar base.
