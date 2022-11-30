abstract type AbstractTimeParallelSolver <: AbstractInitialValueSolver end
abstract type AbstractTimeParallelSolution <: AbstractInitialValueSolution end
abstract type AbstractTimeParallelIterate <: AbstractInitialValueSolution end
abstract type AbstractTimeParallelCache <: AbstractInitialValueCache end
abstract type AbstractTimeParallelParameters <: AbstractInitialValueParameters end

abstract type AbstractPararealParameters <: AbstractTimeParallelParameters end
abstract type AbstractTolerance <: AbstractTimeParallelParameters end
abstract type AbstractWeights <: AbstractTimeParallelParameters end
