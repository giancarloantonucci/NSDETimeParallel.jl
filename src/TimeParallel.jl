module TimeParallel

export TimeParallelSolver
export TimeParallelSolution
export Parareal
export 𝜑₁, 𝜑₂

using Reexport
using Distributed
using ArrowMacros
using LinearAlgebra
@reexport using NSDEBase
using RecipesBase

include("solver.jl")
include("iterate.jl")
include("solution.jl")
include("error.jl")
include("parareal.jl")
include("solve.jl")
include("misc.jl")
include("plot.jl")

end
