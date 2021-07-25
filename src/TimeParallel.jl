module TimeParallel

export TimeParallelSolver, TimeParallelSolution, Parareal
export MovingWindowSolver, MovingWindowSolution, MoWA
export 𝜑₁, 𝜑₂

using Reexport
using Distributed
using ArrowMacros
using LinearAlgebra
@reexport using NSDEBase
using RecipesBase

abstract type TimeParallelSolver <: InitialValueSolver end

include("iterate.jl")
include("solution.jl")
include("error.jl")
include("parareal.jl")
# include("moving.jl")
include("solve.jl")
include("misc.jl")
include("recipes.jl")

end
