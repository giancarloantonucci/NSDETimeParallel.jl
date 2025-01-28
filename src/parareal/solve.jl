function (parareal::Parareal)(cache::PararealCache, solution::PararealSolution, problem::AbstractInitialValueProblem)
    @↓ executionmode = parareal
    if executionmode == "SERIAL"
        parareal_serial!(cache, solution, problem, parareal)
    elseif nprocs() == 1 && !(executionmode == "SERIAL")
        @warn "nprocs() is 1. Default to SERIAL executionmode instead."
        parareal_serial!(cache, solution, problem, parareal)
    elseif nprocs() > 1 && executionmode == "DISTRIBUTED"
        parareal_distributed!(cache, solution, problem, parareal)
    else
        parareal_mpi!(cache, solution, problem, parareal)
    end
    return solution
end

function (parareal::Parareal)(solution::PararealSolution, problem::AbstractInitialValueProblem)
    cache = PararealCache(problem, parareal)
    parareal(cache, solution, problem)
    return solution
end

function (parareal::Parareal)(problem::AbstractInitialValueProblem)
    solution = PararealSolution(problem, parareal)
    parareal(solution, problem)
    return solution
end

# function (parareal::Parareal)(cache::PararealCache, problem::AbstractInitialValueProblem)
#     solution = PararealSolution(problem, parareal)
#     parareal(cache, solution, problem)
#     return solution
# end
