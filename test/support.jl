# Shared fixtures and helpers for this package test suite.
# The truth for every accuracy assertion in this suite is the SERIAL FINE
# SOLVE over the whole span — never the analytic solution. Parareal's job is
# to reproduce the fine solver's answer, errors and all.

const problem = Logistic(0.3, (0.0, 1.0))
const finesolver = RK4(; h=1e-3)
const coarsesolver = RungeKutta4(; h=5e-2)
const finetruth = solve(problem, finesolver)

function chunkgrid(problem, N)
    return [(N - n + 1) / N * problem.tspan[1] + (n - 1) / N * problem.tspan[2] for n in 1:(N + 1)]
end

function make_parareal(;
    N=5, K=N, ϵ=1e-12, coarse=coarsesolver, criterion=RelativeUpdate(; atol=ϵ, rtol=ϵ)
)
    return Parareal(finesolver, coarse; parameters=PararealParameters(; N=N, K=K), criterion)
end
