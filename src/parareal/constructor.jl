"""
    Parareal <: AbstractTimeParallelSolver

A composite type for the parareal algorithm.

# Constructors
```julia
Parareal(finesolver, coarsolver, tolerance, parameters, saveiterates)
Parareal(finesolver, coarsolver; tolerance=Tolerance(), parameters=PararealParameters(), saveiterates=false)
```

## Arguments
- `finesolver :: AbstractInitialValueSolver` : fine solver.
- `coarsolver :: AbstractInitialValueSolver` : coarse solver.
- `parameters :: AbstractPararealParameters` : basic parameters.
- `tolerance :: AbstractTolerance` : tolerance and error mechanism.
- `saveiterates :: Bool` : flags when to save all the Parareal iterates.

# Methods

    (parareal::Parareal)(solution::PararealSolution, problem::AbstractInitialValueProblem; mode::String="SERIAL", saveiterates::Bool=false)
    (parareal::Parareal)(problem::AbstractInitialValueProblem; mode::String="SERIAL", saveiterates::Bool=false)

returns the `solution` of a `problem` using `parareal`. `mode` selects the implementations, `saveiterates` flags when to save all iterates into `solution.iterates`.
"""
struct Parareal{finesolver_T<:AbstractInitialValueSolver, coarsolver_T<:AbstractInitialValueSolver, parameters_T<:AbstractPararealParameters, tolerance_T<:AbstractTolerance, saveiterates_T<:Bool} <: AbstractTimeParallelSolver
    finesolver::finesolver_T
    coarsolver::coarsolver_T
    parameters::parameters_T
    tolerance::tolerance_T
    saveiterates::saveiterates_T
    # function Parareal(finesolver::finesolver_T, coarsolver::coarsolver_T, parameters::parameters_T, tolerance::tolerance_T, saveiterates::saveiterates_T) where {finesolver_T<:AbstractInitialValueSolver, coarsolver_T<:AbstractInitialValueSolver, parameters_T<:AbstractPararealParameters, tolerance_T<:AbstractTolerance, saveiterates_T<:Bool}
    #     finesolver2 = deepcopy(finesolver)
    #     coarsolver2 = deepcopy(coarsolver)
    #     return new{finesolver_T, coarsolver_T, parameters_T, tolerance_T, saveiterates_T}(finesolver2, coarsolver2, parameters, tolerance, saveiterates)
    # end
end

Parareal(finesolver::AbstractInitialValueSolver, coarsolver::AbstractInitialValueSolver; parameters::AbstractPararealParameters = PararealParameters(), tolerance::AbstractTolerance = Tolerance(), saveiterates::Bool = false) = Parareal(finesolver, coarsolver, parameters, tolerance, saveiterates)
