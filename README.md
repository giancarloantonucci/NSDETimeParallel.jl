# NSDETimeParallel.jl

A Julia package implementing time-parallel methods.

[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://giancarloantonucci.github.io/NSDETimeParallel.jl/dev) ![Build Status](https://img.shields.io/github/actions/workflow/status/giancarloantonucci/NSDETimeParallel.jl/CI.yml) ![Coverage Status](https://img.shields.io/codecov/c/github/giancarloantonucci/NSDETimeParallel.jl)

## Installation

<!-- This package is a [registered package](https://juliahub.com/ui/Search?q=NSDETimeParallel&type=packages) compatible with Julia v1.10 and above. From the Julia REPL,

```
]add NSDETimeParallel
``` -->

This package is compatible with Julia v1.10 and above. From the Julia REPL,

```
]add https://github.com/giancarloantonucci/NSDETimeParallel.jl
```

Read the [documentation](https://giancarloantonucci.github.io/NSDETimeParallel.jl/dev) for a complete overview of this package.

## Usage

```julia
using NSDERungeKutta, NSDETimeParallel
using Plots, LaTeXStrings

u0 = [2.0, 3.0, -14.0]
problem = Lorenz(u0, 0.0, 10.0)
solver = Parareal(RK4(h = 1e-3), RK4(h = 1e-1), P = 10)
solution = solve(problem, solver)
plot(solution, xlabel = L"t", label = [L"x" L"y" L"z"])
```

![svg](imgs/lorenz.svg)

## Available methods

This package currently supports only `Parareal`.
