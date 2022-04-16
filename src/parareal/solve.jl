function (parareal::Parareal)(cache::PararealCache, solution::PararealSolution, problem::AbstractInitialValueProblem; mode::String="SERIAL")
    if nprocs() == 1 || mode == "SERIAL"
        parareal_serial!(cache, solution, problem, parareal)
    elseif nprocs() > 1 && mode == "DISTRIBUTED"
        parareal_distributed!(cache, solution, problem, parareal)
    elseif nprocs() > 1 && mode == "MPI"
        parareal_mpi!(cache, solution, problem, parareal)
    else
        error("$mode is not available. Choose either SERIAL, DISTRIBUTED, or MPI.")
    end
    return solution
end

function (parareal::Parareal)(cache::PararealCache, problem::AbstractInitialValueProblem; mode::String="SERIAL", saveiterates::Bool=false)
    solution = PararealSolution(problem, parareal; saveiterates=saveiterates)
    parareal(cache, solution, problem; mode=mode)
    return solution
end

function (parareal::Parareal)(solution::PararealSolution, problem::AbstractInitialValueProblem; mode::String="SERIAL")
    cache = coarseguess(problem, parareal)
    parareal(cache, solution, problem; mode=mode)
    return solution
end

function (parareal::Parareal)(problem::AbstractInitialValueProblem; mode::String="SERIAL", saveiterates::Bool=false)
    solution = PararealSolution(problem, parareal; saveiterates=saveiterates)
    parareal(solution, problem; mode=mode)
    return solution
end
