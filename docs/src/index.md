# TimeParallel.jl

This is the documentation of [TimeParallel.jl](https://github.com/giancarloantonucci/TimeParallel.jl), a Julia package implementing time-parallel methods.

## Manual

```@contents
Depth = 3
```

## API

All exported types and functions are considered part of the public API and thus documented in this manual.

### Abstract types

```@docs
AbstractTimeParallelSolver
AbstractTimeParallelSolution
AbstractTimeParallelIterate
AbstractTimeParallelCache
AbstractTimeParallelParameters
```

### Composite types

```@docs
Parareal
Tolerance
Weights
PararealIterate
PararealSolution
PararealCache
```

### Functions

```@docs
solve
solve!
coarseguess!
```

### Utilities

```@docs
ψ₁
ψ₂
update!
getindex
lastindex
length
resize!
setindex!
```

## Index

```@index
```
