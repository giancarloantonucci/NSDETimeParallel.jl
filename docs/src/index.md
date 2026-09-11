# NSDETimeParallel.jl

This is the documentation of [NSDETimeParallel.jl](https://github.com/giancarloantonucci/NSDETimeParallel.jl), a Julia package implementing time-parallel methods — the Parareal algorithm over one shared skeleton with four execution backends.

## Installation

From the Julia REPL,

```
]add https://github.com/giancarloantonucci/NSDETimeParallel.jl
```

## Getting started

```julia
using NSDERungeKutta, NSDETimeParallel

problem = Lorenz([2.0, 3.0, -14.0], (0.0, 10.0))
parareal = Parareal(RK4(h = 1e-3), RK4(h = 1e-1);       # fine, coarse
                    parameters = PararealParameters(N = 10),
                    tolerance = Tolerance(ϵ = 1e-8))
solution = solve(problem, parareal)
solution.errors        # error per iteration
solution(5.0)          # interpolate
```

The convergence measure lives in `Tolerance`: the unweighted `ψ₁` (default) or the weighted `ψ₂` with [`Weights`](@ref) — the latter is what the zooming strategy of NSDEMovingWindow.jl acts on.

- The [Backends](backends.md) page covers the four execution modes and how to judge speed-ups.
- The [API](api.md) holds the full reference.
