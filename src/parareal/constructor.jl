"""
    Parareal <: AbstractTimeParallelSolver

A composite type for the Parareal algorithm.

# Constructors
```julia
Parareal(finesolver, coarsesolver, parameters, tolerance)
Parareal(finesolver, coarsesolver; parameters=PararealParameters(), tolerance=Tolerance())
```

## Arguments
- `finesolver :: AbstractInitialValueSolver` : fine solver (accurate but expensive).
- `coarsesolver :: AbstractInitialValueSolver` : coarse solver (rough but quick).
- `parameters :: AbstractPararealParameters` : parameters for the correction step.
- `tolerance :: AbstractTolerance` : tolerance and error mechanism.

# Methods

    (parareal::Parareal)(solution::PararealSolution, problem::AbstractInitialValueProblem)
    (parareal::Parareal)(problem::AbstractInitialValueProblem)

returns the `solution` of a `problem` using `parareal`.
"""
struct Parareal{
            finesolver_T<:AbstractInitialValueSolver,
            coarsolver_T<:AbstractInitialValueSolver,
            parameters_T<:AbstractPararealParameters,
            tolerance_T<:AbstractTolerance,
        } <: AbstractTimeParallelSolver
    finesolver::finesolver_T
    coarsesolver::coarsolver_T
    parameters::parameters_T
    tolerance::tolerance_T
end

Parareal(
    finesolver::AbstractInitialValueSolver,
    coarsesolver::AbstractInitialValueSolver;
    parameters::AbstractPararealParameters=PararealParameters(),
    tolerance::AbstractTolerance=Tolerance(),
) = Parareal(finesolver, coarsesolver, parameters, tolerance)
