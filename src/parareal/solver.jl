struct Parareal{objective_T, mode_T, ℱ_T, 𝒢_T, P_T, K_T} <: TimeParallelSolver
    objective::objective_T
    mode::mode_T
    ℱ::ℱ_T
    𝒢::𝒢_T
    P::P_T
    K::K_T
end

function Parareal(ℱ::Function, 𝒢::Function; 𝜑 = 𝜑₁, ϵ = 1e-12, P = 10, K = P, mode = "SERIAL")
    objective = ErrorFunction(𝜑, ϵ)
    return Parareal(objective, mode, ℱ, 𝒢, P, K)
end

function Parareal(finesolver::InitialValueSolver, coarsolver::InitialValueSolver; 𝜑 = 𝜑₁, ϵ = 1e-12, P = 10, K = P, mode = "SERIAL")
    objective = ErrorFunction(𝜑, ϵ)
    function ℱ(problem, u0, tspan)
        subproblem = IVP(problem.rhs, u0, tspan)
        solve(subproblem, finesolver)
    end
    ℱ(problem, u0, t0, tN) = ℱ(problem, u0, (t0, tN))
    function 𝒢(problem, u0, tspan)
        subproblem = IVP(problem.rhs, u0, tspan)
        solve(subproblem, coarsolver)
    end
    𝒢(problem, u0, t0, tN) = 𝒢(problem, u0, (t0, tN))
    return Parareal(objective, mode, ℱ, 𝒢, P, K)
end
