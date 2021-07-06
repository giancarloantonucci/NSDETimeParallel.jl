# TimeParallel

[![Build Status](https://github.com/antonuccig/TimeParallel.jl/workflows/CI/badge.svg)](https://github.com/antonuccig/TimeParallel.jl/actions)
[![Coverage](https://codecov.io/gh/antonuccig/TimeParallel.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/antonuccig/TimeParallel.jl)

## Installation

`TimeParallel` is compatible with Julia `v1.0` and above. From the Julia REPL,
```julia
]add https://github.com/antonuccig/TimeParallel.jl
```

## Usage

```julia
using RungeKutta
using TimeParallel

u0 = [2.0, 3.0, -14.0]
tspan = (0.0, 1.0)
problem = Lorenz(u0, tspan)

fine = Midpoint(h = 1e-3)
coarse = Midpoint(h = 1e-2)
solver = Parareal(fine, coarse, P = 10)

solution = solve(problem, solver)
```

```julia
using Plots
plot(solution)
# savefig("lorenz.svg")
```

![svg](images/lorenz.svg)

## Available methods

`TimeParallel` currently supports only `Parareal`.
