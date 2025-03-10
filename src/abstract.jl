"An abstract type for time-parallel solvers of [`NSDEBase.AbstractInitialValueProblem`](@extref)s."
abstract type AbstractTimeParallelSolver <: AbstractInitialValueSolver end

"An abstract type for time-parallel solutions of [`NSDEBase.AbstractInitialValueProblem`](@extref)s."
abstract type AbstractTimeParallelSolution <: AbstractInitialValueSolution end

"An abstract type for iterates in [`AbstractTimeParallelSolution`](@ref)s."
abstract type AbstractTimeParallelIterate <: AbstractInitialValueSolution end

"An abstract type for caching intermediate computations in [`AbstractTimeParallelSolver`](@ref)s."
abstract type AbstractTimeParallelCache <: AbstractInitialValueCache end

"An abstract type for parameters in [`AbstractTimeParallelSolver`](@ref)s."
abstract type AbstractTimeParallelParameters <: AbstractInitialValueParameters end

"An abstract type for parameters of [`Parareal`](@ref)s."
abstract type AbstractPararealParameters <: AbstractTimeParallelParameters end

"An abstract type for tolerance parameters in [`AbstractTimeParallelSolver`](@ref)s."
abstract type AbstractTolerance <: AbstractTimeParallelParameters end

"An abstract type for proximity-function weights used by [`AbstractTimeParallelSolver`](@ref)s."
abstract type AbstractWeights <: AbstractTimeParallelParameters end
