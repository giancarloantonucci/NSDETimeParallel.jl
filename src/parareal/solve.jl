function (parareal::Parareal)(cache::PararealCache, solution::PararealSolution, problem::AbstractInitialValueProblem, executionmode::String="MPI")
    if executionmode == "SERIAL"
        parareal_serial!(cache, solution, problem, parareal)
    elseif nprocs() == 1 && !(executionmode == "SERIAL")
        @warn "nprocs() is 1. Default to SERIAL executionmode instead."
        parareal_serial!(cache, solution, problem, parareal)
    elseif nprocs() > 1 && executionmode == "DISTRIBUTED"
        parareal_distributed!(cache, solution, problem, parareal)
    elseif nprocs() > 1 && executionmode == "MPI"
        parareal_mpi!(cache, solution, problem, parareal)
    end
    return solution
end

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

# function (parareal::Parareal)(cache::PararealCache, problem::AbstractInitialValueProblem; kwargs...)
#     solution = PararealSolution(problem, parareal)
#     parareal(cache, solution, problem; kwargs...)
#     return solution
# end
