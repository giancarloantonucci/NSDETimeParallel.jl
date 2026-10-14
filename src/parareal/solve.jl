"""
    NSDEBase.initialize_cache(problem, parareal::Parareal) :: PararealCache

returns a reusable [`PararealCache`](@ref) for `problem`, honouring the same
contract NSDERungeKutta honours for its solvers. Build once, then call
`parareal(cache, solution, problem; ...)` repeatedly — e.g. in timing loops,
to report setup costs separately from repeated solve costs. Reuse is not a
zero-allocation guarantee.
"""
function NSDEBase.initialize_cache(problem::AbstractInitialValueProblem, parareal::Parareal)
    return PararealCache(problem, parareal)
end

"""
    NSDEBase.initialize_solution(problem, parareal::Parareal; history=:none) :: PararealSolution

returns a [`PararealSolution`](@ref) for `problem`. Reusable across solves.
"""
function NSDEBase.initialize_solution(
    problem::AbstractInitialValueProblem, parareal::Parareal; history::Symbol=:none
)
    return PararealSolution(problem, parareal; history=history)
end

function history_options(history)
    history in (:none, :interfaces, :iterates) ||
        throw(ArgumentError("history must be :none, :interfaces, or :iterates"))
    return (history === :iterates, history === :interfaces)
end
function (parareal::Parareal)(
    cache::PararealCache,
    solution::PararealSolution,
    problem::AbstractInitialValueProblem;
    backend::AbstractPararealBackend=SerialBackend(),
    history::Symbol=:none,
    directory::Union{Nothing,String}=nothing,
    nocollect::Bool=false,
    initial_weight=nothing,
)
    prepare_problem!(cache, problem, parareal)
    saveiterates, record_interfaces = history_options(history)
    if saveiterates && solution.iterates === nothing
        throw(ArgumentError("construct the solution with history=:iterates"))
    end
    if backend isa PipelinedMPIBackend
        record_interfaces &&
            throw(ArgumentError("interface histories require a bulk-synchronous backend"))
        initial_weight !== nothing &&
            throw(ArgumentError("configure a fixed pipeline weight on the criterion"))
        prepare_chunks!(cache, parareal)
        return parareal!(
            cache, solution, problem, parareal, backend; saveiterates, directory, nocollect
        )
    end
    return parareal!(
        cache,
        solution,
        problem,
        parareal,
        backend;
        saveiterates,
        directory,
        nocollect,
        record_interfaces,
        initial_weight,
    )
end
function (parareal::Parareal)(
    solution::PararealSolution, problem::AbstractInitialValueProblem; kwargs...
)
    return parareal(PararealCache(problem, parareal), solution, problem; kwargs...)
end
function (parareal::Parareal)(
    problem::AbstractInitialValueProblem; history::Symbol=:none, kwargs...
)
    solution = PararealSolution(problem, parareal; history)
    return parareal(PararealCache(problem, parareal), solution, problem; history, kwargs...)
end
function (parareal::Parareal)(
    cache::PararealCache, problem::AbstractInitialValueProblem; history::Symbol=:none, kwargs...
)
    solution = PararealSolution(problem, parareal; history)
    return parareal(cache, solution, problem; history, kwargs...)
end
