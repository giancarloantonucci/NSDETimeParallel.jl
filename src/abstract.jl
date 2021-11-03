"An abstract type for time-parallel solvers."
abstract type AbstractTimeParallelSolver <: AbstractInitialValueSolver end

"An abstract type for the solution of an [`AbstractInitialValueProblem`](@ref) obtained with an [`AbstractTimeParallelSolver`](@ref)."
abstract type AbstractTimeParallelSolution <: AbstractInitialValueSolution end

"An abstract type for the iterate of an [`AbstractTimeParallelSolution`](@ref)."
abstract type AbstractTimeParallelIterate <: AbstractNSDEObject end

"An abstract type for the cache of an [`AbstractTimeParallelSolver`](@ref)."
abstract type AbstractTimeParallelCache <: AbstractNSDEObject end

"An abstract type for the parameters of an [`AbstractTimeParallelSolver`](@ref)."
abstract type AbstractTimeParallelParameters <: AbstractNSDEObject end
