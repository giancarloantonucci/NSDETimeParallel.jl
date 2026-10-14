using NSDETimeParallel
using NSDERungeKutta
using Distributed
using LinearAlgebra
using Serialization
using Test
using Aqua
using RecipesBase: RecipesBase
using NSDEBase: NSDEBase

include("support.jl")

@testset "NSDETimeParallel" begin
    @testset "Aqua" begin
        Aqua.test_all(
            NSDETimeParallel; piracies=(; treat_as_own=[NSDETimeParallel.RecipesBase.recipetype])
        )
    end

    include("cost_models.jl")
    include("criteria.jl")
    include("weights.jl")
    include("parareal.jl")
    include("solutions.jl")
    include("history_and_reuse.jl")
    include("backends.jl")
    include("chunk_schedules.jl")
    include("diagnostics.jl")
    include("problem_reuse.jl")
    include("failures_and_plots.jl")
end
