module TimeParallel

export TimeParallelSolver, TimeParallelSolution, Parareal
export MovingWindowSolver, MovingWindowSolution, MoWA
export 𝜑₁, 𝜑₂

using Reexport
using ArrowMacros
using LinearAlgebra
@reexport using NSDEBase
using RecipesBase

abstract type TimeParallelSolver <: InitialValueSolver end

include("solution.jl")
include("objective.jl")
include("parareal/solver.jl")
include("parareal/solve.jl")
# include("moving.jl")
include("solve.jl")
include("recipes.jl")

end
