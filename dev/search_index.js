var documenterSearchIndex = {"docs":
[{"location":"#NSDETimeParallel.jl","page":"Home","title":"NSDETimeParallel.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This is the documentation of NSDETimeParallel.jl, a Julia package implementing time-parallel methods.","category":"page"},{"location":"#API","page":"Home","title":"API","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"All exported types and functions are considered part of the public API, and thus documented in this manual.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [NSDETimeParallel]","category":"page"},{"location":"#NSDETimeParallel.AbstractPararealParameters","page":"Home","title":"NSDETimeParallel.AbstractPararealParameters","text":"An abstract type for parameters of Parareals.\n\n\n\n\n\n","category":"type"},{"location":"#NSDETimeParallel.AbstractTimeParallelCache","page":"Home","title":"NSDETimeParallel.AbstractTimeParallelCache","text":"An abstract type for caching intermediate computations in AbstractTimeParallelSolvers.\n\n\n\n\n\n","category":"type"},{"location":"#NSDETimeParallel.AbstractTimeParallelIterate","page":"Home","title":"NSDETimeParallel.AbstractTimeParallelIterate","text":"An abstract type for iterates in AbstractTimeParallelSolutions.\n\n\n\n\n\n","category":"type"},{"location":"#NSDETimeParallel.AbstractTimeParallelParameters","page":"Home","title":"NSDETimeParallel.AbstractTimeParallelParameters","text":"An abstract type for parameters in AbstractTimeParallelSolvers.\n\n\n\n\n\n","category":"type"},{"location":"#NSDETimeParallel.AbstractTimeParallelSolution","page":"Home","title":"NSDETimeParallel.AbstractTimeParallelSolution","text":"An abstract type for time-parallel solutions of NSDEBase.AbstractInitialValueProblems.\n\n\n\n\n\n","category":"type"},{"location":"#NSDETimeParallel.AbstractTimeParallelSolver","page":"Home","title":"NSDETimeParallel.AbstractTimeParallelSolver","text":"An abstract type for time-parallel solvers of NSDEBase.AbstractInitialValueProblems.\n\n\n\n\n\n","category":"type"},{"location":"#NSDETimeParallel.AbstractTolerance","page":"Home","title":"NSDETimeParallel.AbstractTolerance","text":"An abstract type for tolerance parameters in AbstractTimeParallelSolvers.\n\n\n\n\n\n","category":"type"},{"location":"#NSDETimeParallel.AbstractWeights","page":"Home","title":"NSDETimeParallel.AbstractWeights","text":"An abstract type for proximity-function weights used by AbstractTimeParallelSolvers.\n\n\n\n\n\n","category":"type"},{"location":"#NSDETimeParallel.Parareal","page":"Home","title":"NSDETimeParallel.Parareal","text":"Parareal <: AbstractTimeParallelSolver\n\nA composite type for the Parareal algorithm.\n\nConstructors\n\nParareal(finesolver, coarsesolver, parameters, tolerance)\nParareal(finesolver, coarsesolver; parameters=PararealParameters(), tolerance=Tolerance())\n\nArguments\n\nfinesolver :: AbstractInitialValueSolver : fine solver (accurate but expensive).\ncoarsesolver :: AbstractInitialValueSolver : coarse solver (rough but quick).\nparameters :: AbstractPararealParameters : parameters for the correction step.\ntolerance :: AbstractTolerance : tolerance and error mechanism.\n\nMethods\n\n(parareal::Parareal)(solution::PararealSolution, problem::AbstractInitialValueProblem)\n(parareal::Parareal)(problem::AbstractInitialValueProblem)\n\nreturns the solution of a problem using parareal.\n\n\n\n\n\n","category":"type"},{"location":"#NSDETimeParallel.PararealIterate","page":"Home","title":"NSDETimeParallel.PararealIterate","text":"PararealIterate <: AbstractTimeParallelIterate\n\nA composite type for a single iterate in a PararealSolution.\n\nConstructors\n\nPararealIterate(chunks::AbstractVector{𝕊}) where 𝕊<:AbstractInitialValueSolution\nPararealIterate(problem::AbstractInitialValueProblem, parareal::Parareal)\n\nFunctions\n\nfirstindex : first index.\ngetindex : get chunk.\nlastindex : last index.\nlength : number of chunks.\nnumchunks : number of chunks.\nsetindex! : set chunk.\n\nMethods\n\n(iterate::PararealIterate)(t::Real)\n\nreturns the value of iterate at t via interpolation.\n\n\n\n\n\n","category":"type"},{"location":"#NSDETimeParallel.PararealParameters","page":"Home","title":"NSDETimeParallel.PararealParameters","text":"PararealParameters <: AbstractPararealParameters\n\nA composite type for the basic parameters of Parareal.\n\nConstructors\n\nPararealParameters(N, K)\nPararealParameters(; N=10, K=N)\n\nArguments\n\nN :: Integer : number of time chunks/processors.\nK :: Integer : maximum number of iterations.\n\n\n\n\n\n","category":"type"},{"location":"#NSDETimeParallel.PararealSolution","page":"Home","title":"NSDETimeParallel.PararealSolution","text":"PararealSolution <: AbstractTimeParallelSolution\n\nA composite type for an AbstractTimeParallelSolution obtained using Parareal.\n\nConstructors\n\nPararealSolution(lastiterate, errors)\nPararealSolution(problem::AbstractInitialValueProblem, parareal::Parareal)\n\nArguments\n\nlastiterate :: PararealIterate\nerrors :: AbstractVector{ℝ} where ℝ<:Real : iteration errors.\n\nFunctions\n\nfirstindex : first index.\ngetindex : get iterate.\nlastindex : last index.\nnumiterates : number of iterates.\nnumchunks : number of chunks of last iterate.\nsetindex! : set iterate.\n\nMethods\n\n(solution::PararealSolution)(t::Real)\n\nreturns the value of solution at t via interpolation.\n\n\n\n\n\n","category":"type"},{"location":"#NSDETimeParallel.Tolerance","page":"Home","title":"NSDETimeParallel.Tolerance","text":"Tolerance <: AbstractTolerance\n\nA composite type for the tolerance mechanism of an time-parallel solver.\n\nConstructors\n\nTolerance(ϵ, ψ, weights)\nTolerance(; ϵ=1e-12, ψ=ψ₁, weights=Weights())\n\nArguments\n\nϵ :: Real : tolerance.\nψ :: Function : error function.\nweights :: Weights : weights for ψ.\n\n\n\n\n\n","category":"type"},{"location":"#NSDETimeParallel.Weights","page":"Home","title":"NSDETimeParallel.Weights","text":"Weights <: AbstractWeights\n\nA composite type for the weights of Tolerance.\n\nConstructors\n\nWeights(; w=1.0, updatew=false)\n\nArguments\n\nw :: Union{AbstractVector{ℝ}, ℝ} where ℝ<:Real : weighting factor for ψ.\nupdatew :: Bool : flags when to update! w using (an approximation of) the Lipschitz function of the fine solver.\n\nFunctions\n\nupdate! : updates w using (an approximation of) the Lipschitz function of the fine solver.\n\n\n\n\n\n","category":"type"},{"location":"#Base.firstindex-Tuple{PararealIterate}","page":"Home","title":"Base.firstindex","text":"firstindex(iterate::PararealIterate)\n\nreturns the first index of iterate.\n\n\n\n\n\n","category":"method"},{"location":"#Base.firstindex-Tuple{PararealSolution}","page":"Home","title":"Base.firstindex","text":"firstindex(solution::PararealSolution)\n\nreturns the first index of solution.\n\n\n\n\n\n","category":"method"},{"location":"#Base.getindex-Tuple{PararealIterate, Integer}","page":"Home","title":"Base.getindex","text":"getindex(iterate::PararealIterate, n::Integer)\n\nreturns the n-th chunk of iterate.\n\n\n\n\n\n","category":"method"},{"location":"#Base.getindex-Tuple{PararealSolution, Integer}","page":"Home","title":"Base.getindex","text":"getindex(solution::PararealSolution, n::Integer)\n\nreturns the n-th chunk of the last iteration of a PararealSolution.\n\n\n\n\n\n","category":"method"},{"location":"#Base.lastindex-Tuple{PararealIterate}","page":"Home","title":"Base.lastindex","text":"lastindex(iterate::PararealIterate)\n\nreturns the last index of iterate.\n\n\n\n\n\n","category":"method"},{"location":"#Base.lastindex-Tuple{PararealSolution}","page":"Home","title":"Base.lastindex","text":"lastindex(solution::PararealSolution)\n\nreturns the last index of solution.\n\n\n\n\n\n","category":"method"},{"location":"#Base.length-Tuple{PararealIterate}","page":"Home","title":"Base.length","text":"length(iterate::PararealIterate)\n\nreturns the number of chunks of iterate.\n\n\n\n\n\n","category":"method"},{"location":"#Base.length-Tuple{PararealSolution}","page":"Home","title":"Base.length","text":"length(solution::PararealSolution)\n\nreturns the number of chunks of solution.\n\n\n\n\n\n","category":"method"},{"location":"#Base.setindex!-Tuple{PararealIterate, AbstractInitialValueSolution, Integer}","page":"Home","title":"Base.setindex!","text":"setindex!(iterate::PararealIterate, value::AbstractInitialValueSolution, n::Integer)\n\nstores value into the n-th chunk of iterate.\n\n\n\n\n\n","category":"method"},{"location":"#Base.setindex!-Tuple{PararealSolution, AbstractInitialValueSolution, Integer}","page":"Home","title":"Base.setindex!","text":"setindex!(solution::PararealSolution, chunk::AbstractInitialValueSolution, n::Integer)\n\nstores an AbstractInitialValueSolution as the n-th chunk of the last iteration of a PararealSolution.\n\n\n\n\n\n","category":"method"},{"location":"#NSDEBase.solve!-Tuple{AbstractTimeParallelSolution, AbstractInitialValueProblem, AbstractTimeParallelSolver}","page":"Home","title":"NSDEBase.solve!","text":"solve!(solution::AbstractTimeParallelSolution, problem, solver::AbstractTimeParallelSolver; kwargs...) :: AbstractTimeParallelSolution\n\ncomputes the solution of problem using solver.\n\n\n\n\n\n","category":"method"},{"location":"#NSDEBase.solve-Tuple{AbstractInitialValueProblem, AbstractTimeParallelSolver}","page":"Home","title":"NSDEBase.solve","text":"solve(problem, solver::AbstractTimeParallelSolver; kwargs...) :: AbstractTimeParallelSolution\n\ncomputes the solution of problem using solver.\n\n\n\n\n\n","category":"method"},{"location":"#NSDETimeParallel.TimeParallelSolution-Tuple{AbstractInitialValueProblem, Parareal}","page":"Home","title":"NSDETimeParallel.TimeParallelSolution","text":"TimeParallelSolution(problem::AbstractInitialValueProblem, parareal::Parareal)\n\nreturns a PararealSolution constructor for the solution of problem with parareal.\n\n\n\n\n\n","category":"method"},{"location":"#NSDETimeParallel.numchunks-Tuple{PararealIterate}","page":"Home","title":"NSDETimeParallel.numchunks","text":"numchunks(iterate::PararealIterate)\n\nreturns the number of chunks of iterate.\n\n\n\n\n\n","category":"method"},{"location":"#NSDETimeParallel.numchunks-Tuple{PararealSolution}","page":"Home","title":"NSDETimeParallel.numchunks","text":"numchunks(solution::PararealSolution)\n\nreturns the number of chunks of solution.\n\n\n\n\n\n","category":"method"},{"location":"#NSDETimeParallel.numiterates-Tuple{PararealSolution}","page":"Home","title":"NSDETimeParallel.numiterates","text":"numiterates(solution::PararealSolution)\n\nreturns the number of iterates of solution.\n\n\n\n\n\n","category":"method"},{"location":"#NSDETimeParallel.parareal_distributed!-Tuple{NSDETimeParallel.PararealCache, PararealSolution, AbstractInitialValueProblem, Parareal}","page":"Home","title":"NSDETimeParallel.parareal_distributed!","text":"Distributed implementation of Parareal.\n\n\n\n\n\n","category":"method"},{"location":"#NSDETimeParallel.parareal_mpi!-Tuple{NSDETimeParallel.PararealCache, PararealSolution, AbstractInitialValueProblem, Parareal}","page":"Home","title":"NSDETimeParallel.parareal_mpi!","text":"MPI implementation of Parareal.\n\n\n\n\n\n","category":"method"},{"location":"#NSDETimeParallel.update!-Union{Tuple{𝕍}, Tuple{ℂ}, Tuple{Weights, AbstractVector{𝕍}, AbstractVector{𝕍}}} where {ℂ<:Number, 𝕍<:AbstractVector{ℂ}}","page":"Home","title":"NSDETimeParallel.update!","text":"update!(weights::Weights, U<:AbstractVector{𝕍}, F<:AbstractVector{𝕍}) where 𝕍<:AbstractVector{ℂ} where ℂ<:Number\n\nupdates weights.w based on U and F.\n\n\n\n\n\n","category":"method"},{"location":"#NSDETimeParallel.ψ₁-Tuple{Any, Any, Any}","page":"Home","title":"NSDETimeParallel.ψ₁","text":"standard error function.\n\n\n\n\n\n","category":"method"},{"location":"#NSDETimeParallel.ψ₂-Tuple{Any, Any, Any}","page":"Home","title":"NSDETimeParallel.ψ₂","text":"weighted error function.\n\n\n\n\n\n","category":"method"}]
}
