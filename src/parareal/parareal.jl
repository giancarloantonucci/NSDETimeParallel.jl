"""
    Parareal <: AbstractTimeParallelSolver

A composite type for the parareal algorithm.

# Constructors
```julia
Parareal(finesolver, coarsolver, P, K, control)
Parareal(finesolver, coarsolver; P=10, K=P, control=ErrorControl())
```

# Arguments
- `finesolver :: AbstractInitialValueSolver` : fine solver.
- `coarsolver :: AbstractInitialValueSolver` : coarse solver.
- `P :: Integer` : number of time chunks.
- `K :: Integer` : maximum number of iterations.
- `control :: ErrorControl` : tolerance and error mechanism.

# Functions
- [`show`](@ref) : shows name and contents.
- [`summary`](@ref) : shows name.
"""
struct Parareal{finesolver_T, coarsolver_T, P_T, K_T, control_T} <: AbstractTimeParallelSolver
    finesolver::finesolver_T
    coarsolver::coarsolver_T
    P::P_T
    K::K_T
    control::control_T
end

function Parareal(finesolver, coarsolver; P=10, K=P, control=ErrorControl())
    return Parareal(finesolver, coarsolver, P, K, control)
end
