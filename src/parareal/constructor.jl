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
- `saveiterates :: Bool` : save all the Parareal iterates if `true`.

# Methods

    (parareal::Parareal)(solution::PararealSolution, problem::AbstractInitialValueProblem; mode::String="SERIAL", saveiterates::Bool=false)
    (parareal::Parareal)(problem::AbstractInitialValueProblem; mode::String="SERIAL", saveiterates::Bool=false)

returns the `solution` of a `problem` using `parareal`. `mode` selects the implementation, `saveiterates` flags whether to save all iterates in `solution.iterates`.
"""
struct Parareal{finesolver_T<:AbstractInitialValueSolver, coarsolver_T<:AbstractInitialValueSolver, parameters_T<:AbstractPararealParameters, tolerance_T<:AbstractTolerance, saveiterates_T<:Bool} <: AbstractTimeParallelSolver
    finesolver::finesolver_T
    coarsesolver::coarsolver_T
    parameters::parameters_T
    tolerance::tolerance_T
    saveiterates::saveiterates_T
    # function Parareal(finesolver::finesolver_T, coarsesolver::coarsolver_T, parameters::parameters_T, tolerance::tolerance_T, saveiterates::saveiterates_T) where {finesolver_T<:AbstractInitialValueSolver, coarsolver_T<:AbstractInitialValueSolver, parameters_T<:AbstractPararealParameters, tolerance_T<:AbstractTolerance, saveiterates_T<:Bool}
    #     finesolver2 = deepcopy(finesolver)
    #     coarsolver2 = deepcopy(coarsesolver)
    #     return new{finesolver_T, coarsolver_T, parameters_T, tolerance_T, saveiterates_T}(finesolver2, coarsolver2, parameters, tolerance, saveiterates)
    # end
end

Parareal(finesolver::AbstractInitialValueSolver, coarsesolver::AbstractInitialValueSolver; parameters::AbstractPararealParameters = PararealParameters(), tolerance::AbstractTolerance = Tolerance(), saveiterates::Bool = false) = Parareal(finesolver, coarsesolver, parameters, tolerance, saveiterates)
