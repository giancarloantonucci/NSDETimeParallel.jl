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
import RungeKutta
using RecipesBase

include("solver.jl")
include("iterate.jl")
include("solution.jl")
include("error.jl")
include("parareal.jl")
include("solve.jl")
include("utils.jl")
include("plot.jl")

end
