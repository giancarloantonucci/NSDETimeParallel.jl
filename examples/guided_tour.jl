# Executable companion to docs/src/guide.md.
# Run in an environment containing this package and the dependencies used below.

using NSDERungeKutta, NSDETimeParallel, LinearAlgebra
problem = Logistic(0.3, (0.0, 2.0))
fine, coarse = RK4(; h=0.005), RK4(; h=0.05)
par = Parareal(fine, coarse; N=4, criterion=RelativeUpdate(; atol=1e-12, rtol=1e-9))
r = solve(problem, par; history=:interfaces)
@assert r.diagnostics.status in (:converged, :exact)
@assert r.diagnostics.fine_batches == r.diagnostics.sweeps + 1
@assert r.errors === r.diagnostics.residuals
r.diagnostics

full = remake(par; criterion=NoStop())
e = solve(problem, full; history=:interfaces)
Ustar = reference(problem, full)
@assert e.diagnostics.status == :exact
@assert e.diagnostics.fine_batches == 4
@assert e.diagnostics.sweeps == 3
@assert boundaryvalues(e.lastiterate) == Ustar
for (k, U) in enumerate(e.diagnostics.interfaces)
    @assert U[1:k] == Ustar[1:k]
end
nothing

chaos = Lorenz([2.0, 3.0, -14.0], (0.0, 2.0))
F = RK4(; h=0.002)
chunk = remake(chaos; tspan=(0.0, 0.5))
fine_end = solve(chunk, F).u[end]
for hG in (0.05, 0.02, 0.01)
    coarse_end = solve(chunk, RK4(; h=hG)).u[end]
    println((hG=hG, chunk_error=norm(coarse_end-fine_end)))
end
p = Parareal(F, RK4(; h=0.02); N=4, criterion=RelativeUpdate(; rtol=1e-8))
a = solve(chaos, p; backend=ThreadsBackend())
ref = reference(chaos, p)
interface_error = maximum(norm.(boundaryvalues(a.lastiterate) .- ref))
ζ = costratio(chaos, p)
model = theoretical_speedup(a, ζ)
budget = iteration_budget(2.0, 4, ζ)
@assert isfinite(interface_error) && isfinite(model)
@assert 0 <= budget <= 4
println((
    status=a.diagnostics.status,
    batches=a.diagnostics.fine_batches,
    modelled_speedup=model,
    interface_error,
    max_seam=maxseam(a),
))
nothing

weighted = remake(p; criterion=WeightedDefect(; ϵ=1e-8, weight=FixedWeight(exp(0.9))))  # an illustrative rate per unit time
wrun = solve(chaos, weighted)
secant = remake(p; criterion=WeightedDefect(; ϵ=1e-8, weight=SecantWeight()))
srun = solve(chaos, secant)
@assert !isempty(srun.diagnostics.weights)
measured_weight = last(srun.diagnostics.weights)
weighted_error = Wnorm(srun, reference(chaos, secant))
@assert isfinite(weighted_error) && measured_weight >= 1
maximum_rule = WeightedDefect(; ϵ=1e-8, weight=FixedWeight(1), p=Inf)
nothing

saved = solve(problem, par; history=:iterates)
@assert saved.iterates !== nothing
cache = initialize_cache(problem, par)
reused = initialize_solution(problem, par)
new_problem = remake(problem; u0=[0.4])
par(cache, reused, new_problem)
@assert boundaryvalues(reused.lastiterate)[1] == [0.4]
larger = remake(par; N=8)
@assert larger.parameters.K == 8
nothing
