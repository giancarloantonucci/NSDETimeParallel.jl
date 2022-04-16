"""
    Parareal <: AbstractTimeParallelSolver

A composite type for the parareal algorithm.

# Constructors
```julia
Parareal(finesolver, coarsolver, tolerance, P, K)
Parareal(finesolver, coarsolver; tolerance=Tolerance(), P=10, K=P)
```

## Arguments
- `finesolver :: AbstractInitialValueSolver` : fine solver.
- `coarsolver :: AbstractInitialValueSolver` : coarse solver.
- `tolerance :: AbstractTolerance` : tolerance and error mechanism.
- `P :: Integer` : number of time chunks.
- `K :: Integer` : maximum number of iterations.

# Methods

    (parareal::Parareal)(solution::PararealSolution, problem::AbstractInitialValueProblem; mode::String="SERIAL")
    (parareal::Parareal)(problem::AbstractInitialValueProblem; mode::String="SERIAL", saveiterates::Bool=false)

returns the `solution` of a `problem` using `parareal`; `mode` selects the implementations, `saveiterates` is a flag to save all iterates into `solution.alliterates`.
"""
mutable struct Parareal{finesolver_T<:AbstractInitialValueSolver, coarsolver_T<:AbstractInitialValueSolver, tolerance_T<:AbstractTolerance, P_T<:Integer, K_T<:Integer} <: AbstractTimeParallelSolver
    finesolver::finesolver_T
    coarsolver::coarsolver_T
    tolerance::tolerance_T
    P::P_T
    K::K_T
end
Parareal(finesolver::AbstractInitialValueSolver, coarsolver::AbstractInitialValueSolver; tolerance::AbstractTolerance=Tolerance(), P::Integer=10, K::Integer=P) = Parareal(finesolver, coarsolver, tolerance, P, K)
