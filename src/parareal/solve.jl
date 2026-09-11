# src/parareal/solve.jl

"""
    NSDEBase.initialize_cache(problem, parareal::Parareal) :: PararealCache

returns a reusable [`PararealCache`](@ref) for `problem`, honouring the same
contract NSDERungeKutta honours for its solvers. Build once, then call
`parareal(cache, solution, problem; ...)` repeatedly — e.g. in timing loops,
where rebuilding the cache per solve measures the allocator, not the algorithm.
"""
NSDEBase.initialize_cache(problem::AbstractInitialValueProblem, parareal::Parareal) = PararealCache(problem, parareal)

"""
    NSDEBase.initialize_solution(problem, parareal::Parareal; saveiterates=false) :: PararealSolution

returns a [`PararealSolution`](@ref) for `problem`. Reusable across solves.
"""
NSDEBase.initialize_solution(problem::AbstractInitialValueProblem, parareal::Parareal; saveiterates::Bool=false) = PararealSolution(problem, parareal; saveiterates)

_backend_from_mode(mode::AbstractString) =
    uppercase(mode) == "SERIAL"      ? SerialBackend()       :
    uppercase(mode) == "THREADS"     ? ThreadsBackend()      :
    uppercase(mode) == "DISTRIBUTED" ? DistributedBackend()  :
    uppercase(mode) == "MPI"         ? MPIBackend()          :
    uppercase(mode) == "PIPELINED"   ? PipelinedMPIBackend() :
    throw(ArgumentError("Unknown Parareal mode \"$mode\"; use \"SERIAL\", \"THREADS\", \"DISTRIBUTED\", \"MPI\" or \"PIPELINED\", or pass a backend object."))
    # (The old dispatcher SILENTLY returned the unsolved solution on a typo.)

function (parareal::Parareal)(cache::PararealCache, solution::PararealSolution, problem::AbstractInitialValueProblem;
        mode::String="SERIAL", backend::Union{Nothing,AbstractPararealBackend}=nothing,
        saveiterates::Bool=false, directory::Union{Nothing,String}=nothing, nocollect::Bool=false)
    b = backend === nothing ? _backend_from_mode(mode) : backend
    return parareal!(cache, solution, problem, parareal, b; saveiterates, directory, nocollect)
end

function (parareal::Parareal)(solution::PararealSolution, problem::AbstractInitialValueProblem; kwargs...)
    cache = PararealCache(problem, parareal)
    parareal(cache, solution, problem; kwargs...)
    return solution
end

function (parareal::Parareal)(problem::AbstractInitialValueProblem; saveiterates::Bool=false, kwargs...)
    solution = PararealSolution(problem, parareal; saveiterates)
    parareal(solution, problem; saveiterates, kwargs...)
    return solution
end

function (parareal::Parareal)(cache::PararealCache, problem::AbstractInitialValueProblem; saveiterates::Bool=false, kwargs...)
    solution = PararealSolution(problem, parareal; saveiterates)
    parareal(cache, solution, problem; saveiterates, kwargs...)
    return solution
end
