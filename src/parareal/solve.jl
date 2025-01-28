function (parareal::Parareal)(cache::PararealCache, solution::PararealSolution, problem::AbstractInitialValueProblem;
    directory::String="results", mode::String="DISTRIBUTED", saveiterates::Bool=false)
    if mode == "DISTRIBUTED"
        parareal_distributed!(cache, solution, problem, parareal; directory, saveiterates)
    elseif mode == "MPI"
        parareal_mpi!(cache, solution, problem, parareal; directory, saveiterates)
    end
    return solution
end

# function (parareal::Parareal)(cache::PararealCache, solution::PararealSolution, problem::AbstractInitialValueProblem;
#         directory::String="results", saveiterates::Bool=false)
#     parareal_mpi!(cache, solution, problem, parareal; directory, saveiterates)
#     return solution
# end

function (parareal::Parareal)(solution::PararealSolution, problem::AbstractInitialValueProblem; kwargs...)
    cache = PararealCache(problem, parareal)
    parareal(cache, solution, problem; kwargs...)
    return solution
end

function (parareal::Parareal)(problem::AbstractInitialValueProblem; kwargs...)
    solution = PararealSolution(problem, parareal)
    parareal(solution, problem; kwargs...)
    return solution
end

function (parareal::Parareal)(cache::PararealCache, problem::AbstractInitialValueProblem; kwargs...)
    solution = PararealSolution(problem, parareal)
    parareal(cache, solution, problem; kwargs...)
    return solution
end
