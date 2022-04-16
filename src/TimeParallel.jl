module TimeParallel

using Reexport
using Distributed
using ArrowMacros
using LinearAlgebra
@reexport using NSDEBase
using RecipesBase

include("abstract.jl")
include("utils.jl")
include("weights.jl")
include("tolerance.jl")
include("parareal/constructor.jl")
include("parareal/cache.jl")
include("parareal/iterate.jl")
include("parareal/solution.jl")
include("parareal/coarseguess.jl")
include("parareal/serial.jl")
include("parareal/distributed.jl")
include("parareal/mpi.jl")
include("parareal/solve.jl")
include("solve.jl")
include("plot.jl")

export AbstractTimeParallelSolver
export AbstractTimeParallelSolution
export AbstractTimeParallelIterate
export AbstractTimeParallelParameters

export Parareal
export PararealIterate
export PararealSolution

export TimeParallelSolution
export Tolerance, Weights

export coarseguess, coarseguess!
export ψ₁, ψ₂
export fulllength

end
