function (parareal::Parareal)(cache::PararealCache, solution::PararealSolution, problem::AbstractInitialValueProblem;
    directory::String="results", mode::String="DISTRIBUTED", saveiterates::Bool=false, nocollect::Bool=false)
    if mode == "DISTRIBUTED"
        parareal_distributed!(cache, solution, problem, parareal; directory, saveiterates, nocollect)
    elseif mode == "MPI"
        parareal_mpi!(cache, solution, problem, parareal; directory, saveiterates, nocollect)
    end
    return solution
end

function (parareal::Parareal)(solution::PararealSolution, problem::AbstractInitialValueProblem; kwargs...)
    cache = PararealCache(problem, parareal)
    parareal(cache, solution, problem; kwargs...)
    return solution
end

function (parareal::Parareal)(problem::AbstractInitialValueProblem; saveiterates::Bool=false, kwargs...)
    solution = PararealSolution(problem, parareal; saveiterates)
    parareal(solution, problem; kwargs...)
    return solution
end

function (parareal::Parareal)(cache::PararealCache, problem::AbstractInitialValueProblem; saveiterates::Bool=false, kwargs...)
    solution = PararealSolution(problem, parareal; saveiterates)
    parareal(cache, solution, problem; kwargs...)
    return solution
end
