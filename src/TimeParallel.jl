module TimeParallel

using Distributed
using ArrowMacros
using LinearAlgebra
using NSDEBase
using RungeKutta
using RecipesBase

include("abstract.jl")
include("weights.jl")
include("error.jl")
include("parareal/parareal.jl")
include("parareal/cache.jl")
include("parareal/iterate.jl")
include("parareal/solution.jl")
include("parareal/coarse.jl")
include("parareal/serial.jl")
include("parareal/distributed.jl")
include("parareal/mpi.jl")
include("parareal/solve.jl")
include("solve.jl")
include("plotrecipes.jl")

export AbstractTimeParallelSolver
export AbstractTimeParallelSolution
export Parareal
export PararealSolution
export ErrorControl, ErrorWeights
export ψ₁, ψ₂

end
