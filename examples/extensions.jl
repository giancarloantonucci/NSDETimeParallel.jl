# Companion to docs/src/extensions.md.

using NSDERungeKutta, NSDETimeParallel, LinearAlgebra
struct MaxJump <: AbstractCriterion
    tolerance::Float64
end
function NSDETimeParallel.residual(c::MaxJump, cache, k)
    return maximum(norm(cache.U[n] - cache.F[n]) for n in eachindex(cache.U))
end
NSDETimeParallel.threshold(c::MaxJump) = c.tolerance
NSDETimeParallel.accept(c::MaxJump, value) = isfinite(value) && value <= c.tolerance
problem = Logistic(0.3, (0.0, 2.0))
p = Parareal(RK4(; h=0.005), RK4(; h=0.05); N=4, criterion=MaxJump(1e-8))
r = solve(problem, p)
@assert r.diagnostics.status in (:converged, :exact)
nothing

struct Wrapped{S} <: AbstractInitialValueSolver
    inner::S
end
function NSDEBase.initialize_cache(problem::AbstractInitialValueProblem, s::Wrapped)
    return initialize_cache(problem, s.inner)
end
function NSDEBase.initialize_solution(problem::AbstractInitialValueProblem, s::Wrapped; kw...)
    return initialize_solution(problem, s.inner; kw...)
end
function NSDEBase.solve!(cache, sol, problem::AbstractInitialValueProblem, s::Wrapped)
    return solve!(cache, sol, problem, s.inner)
end
function NSDEBase.solve(problem::AbstractInitialValueProblem, s::Wrapped; kw...)
    return solve(problem, s.inner; kw...)
end
function NSDEBase.chunk_solver(s::Wrapped, duration::Real)
    return Wrapped(NSDEBase.chunk_solver(s.inner, duration))
end
pw = Parareal(Wrapped(RK4(; h=0.007)), RK4(; h=0.05); N=4, criterion=NoStop())
rw = solve(problem, pw)
@assert boundaryvalues(rw.lastiterate) == reference(problem, pw)
nothing
