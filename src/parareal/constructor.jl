"""
    Parareal <: AbstractTimeParallelSolver

A composite type for the Parareal algorithm.

# Constructors
```julia
Parareal(finesolver, coarsesolver, tolerance, parameters, saveiterates)
Parareal(finesolver, coarsesolver; tolerance=Tolerance(), parameters=PararealParameters(), saveiterates=false)
```

## Arguments
- `finesolver :: AbstractInitialValueSolver` : fine solver (accurate but expensive).
- `coarsesolver :: AbstractInitialValueSolver` : coarse solver (rough but quick).
- `parameters :: AbstractPararealParameters` : parameters for the correction step.
- `tolerance :: AbstractTolerance` : tolerance and error mechanism.
- `executionmode :: String` : algorithm's kind (`SERIAL`, `DISTRIBUTED`, `MPI`, etc.).
- `saveiterates :: Bool` : save all the Parareal iterates if `true`.

# Methods

    (parareal::Parareal)(solution::PararealSolution, problem::AbstractInitialValueProblem; executionmode::String="SERIAL", saveiterates::Bool=false)
    (parareal::Parareal)(problem::AbstractInitialValueProblem; executionmode::String="SERIAL", saveiterates::Bool=false)

returns the `solution` of a `problem` using `parareal`. `executionmode` selects the implementation, `saveiterates` flags whether to save all iterates in `solution.iterates`.
"""
struct Parareal{
            finesolver_T<:AbstractInitialValueSolver,
            coarsolver_T<:AbstractInitialValueSolver,
            parameters_T<:AbstractPararealParameters,
            tolerance_T<:AbstractTolerance,
            executionmode_T<:String,
            saveiterates_T<:Bool,
        } <: AbstractTimeParallelSolver
    finesolver::finesolver_T
    coarsesolver::coarsolver_T
    parameters::parameters_T
    tolerance::tolerance_T
    executionmode::executionmode_T
    saveiterates::saveiterates_T
end

Parareal(
    finesolver::AbstractInitialValueSolver,
    coarsesolver::AbstractInitialValueSolver;
    parameters::AbstractPararealParameters=PararealParameters(),
    tolerance::AbstractTolerance=Tolerance(),
    executionmode::String="SERIAL",
    saveiterates::Bool=false
) = Parareal(finesolver, coarsesolver, parameters, tolerance, executionmode, saveiterates)
# ) = Parareal(finesolver, coarsesolver, parameters, tolerance, saveiterates)
