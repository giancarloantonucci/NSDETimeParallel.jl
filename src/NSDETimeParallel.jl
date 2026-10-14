module NSDETimeParallel

using Distributed
using MPI
using Serialization

using Reexport
using ArrowMacros
using LinearAlgebra
@reexport using NSDEBase
using RecipesBase

include("abstract.jl")
include("weights.jl")
include("criteria.jl")
include("parareal/parameters.jl")
include("parareal/constructor.jl")
include("parareal/cache.jl")
include("parareal/iterate.jl")
include("parareal/solution.jl")
include("parareal/backends.jl")
include("parareal/parareal.jl")
include("parareal/pipelined.jl")
include("parareal/solve.jl")
include("solve.jl")
include("utils.jl")
include("plots_recipes.jl")

export AbstractTimeParallelSolver
export AbstractTimeParallelSolution
export AbstractTimeParallelIterate
export AbstractTimeParallelParameters

export Parareal
export PararealIterate
export PararealSolution
export PararealParameters
export AbstractCriterion,
    RelativeUpdate, WeightedDefect, MeanRelativeUpdate, NoStop, FrontierRelativeUpdate
export AbstractWeightRule, FixedWeight, SecantWeight, NeighbourWeight
export residual, accept, reference

export PararealDiagnostics
export Wnorm, collect_iterates!, numchunks, numiterates
export boundarytimes, boundaryvalues, flatten
export seams, maxseam
export AbstractPararealBackend,
    SerialBackend, ThreadsBackend, DistributedBackend, MPIBackend, PipelinedMPIBackend
export theoretical_speedup, costratio
export iteration_budget, contractionrate
export shiftwindow!

end
