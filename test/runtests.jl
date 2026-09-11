using NSDETimeParallel
using NSDERungeKutta
using Distributed
using LinearAlgebra
using Serialization
using Test
using Aqua
import RecipesBase
import NSDEBase

# The truth for every accuracy assertion in this suite is the SERIAL FINE
# SOLVE over the whole span — never the analytic solution. Parareal's job is
# to reproduce the fine solver's answer, errors and all.

const problem = Logistic(0.3, (0.0, 1.0))
const finesolver = RK4(h=1e-3)
const coarsesolver = RungeKutta4(h=5e-2)
const finetruth = solve(problem, finesolver)

chunkgrid(problem, N) = [(N - n + 1) / N * problem.tspan[1] + (n - 1) / N * problem.tspan[2] for n = 1:N+1]

make_parareal(; N=5, K=N, ϵ=1e-12, coarse=coarsesolver) =
    Parareal(finesolver, coarse; parameters=PararealParameters(N=N, K=K), tolerance=Tolerance(ϵ=ϵ))

@testset "NSDETimeParallel" begin

@testset "Aqua" begin
    Aqua.test_all(NSDETimeParallel; persistent_tasks=false, piracies=(; treat_as_own=[NSDETimeParallel.RecipesBase.recipetype]))
end

@testset "theoretical_speedup" begin
    @test theoretical_speedup(1, 10, 0.0) ≈ 10.0            # one sweep, free coarse: ideal
    @test theoretical_speedup(2, 10, 0.1) ≈ 10 / (2 * 2 * 0.95) # by hand
    @test theoretical_speedup(2, 10, 0.1) < theoretical_speedup(1, 10, 0.1) # more sweeps, less gain
    @test theoretical_speedup(2, 10, 0.1; wallclock=true) ≈ 10 / (2 * 2) # no work credit
    @test theoretical_speedup(2, 10, 0.1; wallclock=true) < theoretical_speedup(2, 10, 0.1)
    @test theoretical_speedup(2, 2, 0.0; wallclock=true) ≈ 1.0 # K = N cannot win on wall clock
    # The Aubanel (pipelined) ceiling: the coarse chain is hidden behind fine
    # work, so ζN no longer multiplies K.
    @test theoretical_speedup(2, 10, 0.1; estimate=:aubanel) ≈ 10 / (2 * 0.95 * 1.1 + 0.1 * 9) # by hand
    @test theoretical_speedup(4, 100, 0.01; estimate=:aubanel) >
          theoretical_speedup(4, 100, 0.01)                    # pipelining strictly helps at ζ > 0
    @test theoretical_speedup(4, 100, 0.0; estimate=:aubanel) ≈
          theoretical_speedup(4, 100, 0.0)                     # free coarse: nothing to hide
    @test theoretical_speedup(2, 10, 0.1; wallclock=true) ==
          theoretical_speedup(2, 10, 0.1; estimate=:wallclock) # back-compat spelling
    @test_throws ArgumentError theoretical_speedup(2, 10, 0.1; estimate=:naive)
    @test_throws ArgumentError theoretical_speedup(0, 10, 0.1)
end

@testset "PararealParameters validation" begin
    @test_throws ArgumentError PararealParameters(N=0)          # default K = N = 0
    @test_throws ArgumentError PararealParameters(N=5, K=0)
    @test_throws ArgumentError PararealParameters(-1, 3)
    @test_logs (:warn, r"dead work") PararealParameters(N=2, K=5)
    @test PararealParameters(N=5, K=5) isa PararealParameters   # K = N: silent
end

@testset "costratio" begin
    ζ = costratio(problem, make_parareal(N=5))
    @test isfinite(ζ) && ζ > 0
    @test ζ < 1                        # the coarse solver takes 50× fewer steps here
    @test_throws ArgumentError costratio(problem, make_parareal(N=5); repeats=0)
end

@testset "iteration_budget" begin
    # The improved estimate must reproduce tab:sensibleKs verbatim
    # (ζ there is the step-ratio 1/ξ; here it is an INPUT, so that's exact).
    for (N, ζ, S, K) in ((10, 0.1, 2, 2), (10, 0.1, 3, 1), (10, 0.1, 4, 1),
                         (10, 0.01, 2, 6), (10, 0.01, 3, 3), (10, 0.01, 4, 2),
                         (100, 0.1, 2, 4), (100, 0.1, 3, 3), (100, 0.1, 4, 2),
                         (100, 0.01, 2, 29), (100, 0.01, 3, 18), (100, 0.01, 4, 13))
        @test iteration_budget(S, N, ζ) == K
    end
    @test iteration_budget(2, 100, 0.01; estimate=:naive) == 25 # tab anchor, naive form
    # The round-trip identity: the budget is the LARGEST K meeting the
    # target — one more sweep must fall below it.
    for (S, N, ζ) in ((2, 100, 0.01), (3, 64, 0.1), (2, 50, 0.05))
        K = iteration_budget(S, N, ζ)
        @test theoretical_speedup(K, N, ζ) ≥ S > theoretical_speedup(K + 1, N, ζ)
    end
    @test iteration_budget(1e6, 10, 0.1) == 0                 # unreachable target
    # A target BELOW the K = N speed-up is met by every sweep count: the
    # budget is N, the algorithm's own limit. The quadratic's discriminant is
    # negative here, and the old inversion read that as "infeasible" → 0.
    @test theoretical_speedup(10, 10, 0.0) > 0.1
    @test iteration_budget(0.1, 10, 0.0) == 10
    @test iteration_budget(0.1, 10, 0.0; estimate=:naive) == 10 # clamped too, never past N
    @test iteration_budget(1e-300, 10, 0.0; estimate=:naive) == 10 # clamped BEFORE the Int conversion
    @test iteration_budget(1e-300, 10, 0.0) == 10
    @test iteration_budget(1.0, 1, 0.0) == 1                  # N = 1: one sweep, speed-up 1
    # Exhaustive round trip on a grid: the budget is exactly the largest
    # feasible K in 0:N, no more and no less.
    for N in (1, 2, 7, 64, 500), ζ in (0.0, 0.02, 0.3), S in (0.05, 0.9, 1.0, 1.7, 4.0, 60.0)
        K = iteration_budget(S, N, ζ)
        @test 0 ≤ K ≤ N
        K ≥ 1 && @test theoretical_speedup(K, N, ζ) ≥ S
        K < N && @test theoretical_speedup(K + 1, N, ζ) < S
    end
    @test iteration_budget(2, 10, 0.1; estimate=:naive) ==
          floor(Int, 10 / (2 * (1 + 0.1 * 10)))               # naive closed form
    @test iteration_budget(2, 10, 0.1; estimate=:naive) ≤ iteration_budget(2, 10, 0.1) # credit only helps
    @test_throws ArgumentError iteration_budget(0, 10, 0.1)
    @test_throws ArgumentError iteration_budget(2, 10, 0.1; estimate=:pipelined)
end

@testset "contractionrate" begin
    @test contractionrate([1e-1, 1e-3, 1e-5]) ≈ 1e-2          # exact geometric trace
    @test contractionrate([1e-1, 1e-3, 1e-5, 0.0]) ≈ 1e-2     # finite-termination zero excluded
    @test isnan(contractionrate([1e-6]))                       # no observable rate
    @test isnan(contractionrate([0.0, 0.0]))
    solution = solve(problem, make_parareal(N=5))
    β = contractionrate(solution)                              # solution dispatch
    @test isnan(β) || (isfinite(β) && β > 0)
end

@testset "seams and maxseam" begin
    full = solve(problem, make_parareal(N=4))                  # converged run
    @test length(seams(full)) == 3
    @test maxseam(full) == maximum(seams(full))
    @test maxseam(full) < 1e-9                                 # boundaries handed over
    # At ϵ = 0, K = N finite termination makes every START exact, and it is
    # exact in FLOATING POINT too: the newly exact interface is assigned from
    # the fine output, and the correction is evaluated as F + (Gnew − Gold),
    # which returns F bitwise when the coarse values cancel bitwise. (An
    # earlier version computed (Gnew + F) − Gold, lost an ulp, and the test
    # was loosened to 1e-13 with a comment blaming floating point. The
    # ordering was to blame.) The other ingredient is that a fine step which
    # divides the chunk ends ON the chunk boundary, so `solution(T[n+1])` is
    # the final node and not an interpolation — see the chunk-boundary test.
    exact = solve(problem, make_parareal(N=4, K=4, ϵ=0.0))
    @test maxseam(exact) == 0.0
    early = solve(problem, make_parareal(N=4, K=1, ϵ=0.0))     # one sweep: visible seams
    @test all(≥(0), seams(early))
    @test maxseam(early) > maxseam(full)
    single = solve(problem, make_parareal(N=1, K=1))
    @test isempty(seams(single))
    @test maxseam(single) == 0.0
end

@testset "δ is a measurement safety factor, not a ratchet" begin
    # `update!` is internal (deliberately unexported: the name is too generic
    # to claim), so the tests call it qualified.
    frozen = Weights(w=2.0, updatew=false, δ=0.5)
    U = [[0.0], [1.0], [2.0], [3.0]]
    F = [[0.0], [1.0], [3.0], [7.0]]
    T = [0.0, 1.0, 2.0, 3.0]
    NSDETimeParallel.update!(frozen, U, F, T)
    NSDETimeParallel.update!(frozen, U, F, T)
    @test frozen.w == 2.0                  # updatew=false: strict no-op (old code halved it twice)
    measured = Weights(updatew=true, δ=0.5)
    NSDETimeParallel.update!(measured, U, F, T)
    @test measured.w ≈ 8.0                 # max per-time ratio 4, inflated by 1/δ = 2
end

@testset "update! never poisons a valid weight" begin
    T = [0.0, 1.0, 2.0, 3.0]
    # Three coincident starts and coincident fine outputs — an exact
    # equilibrium, or an already-exact prefix — used to give 0/0 = NaN and
    # `max(2.0, NaN) == NaN`, after which ψ₂ was blind for the rest of the run.
    still = Weights(w=2.0, updatew=true)
    U = [[1.0], [1.0], [1.0], [1.0]]
    NSDETimeParallel.update!(still, U, U, T)
    @test still.w == 2.0                   # nothing measurable: previous weight kept
    # A zero denominator with a nonzero numerator (Inf) is skipped as well:
    inf = Weights(w=2.0, updatew=true)
    NSDETimeParallel.update!(inf, [[1.0], [1.0], [1.0], [1.0]], [[0.0], [1.0], [3.0], [7.0]], T)
    @test inf.w == 2.0
    # A usable ratio next to an unusable one is still measured:
    mixed = Weights(w=1.0, updatew=true)
    NSDETimeParallel.update!(mixed, [[0.0], [1.0], [1.0], [3.0]], [[0.0], [1.0], [3.0], [7.0]], T)
    @test mixed.w ≈ 2.0                    # i = 2 gives (3 − 1)/(1 − 0) = 2; i = 3 is 0/0 and skipped
    @test isfinite(mixed.w)
    # Garbage in the constructor is refused rather than stored:
    @test_throws ArgumentError Weights(δ=0.0)
    @test_throws ArgumentError Weights(δ=Inf)
    @test_throws ArgumentError Weights(w=NaN)
    @test_throws ArgumentError Weights(w=0.0)
    @test_throws ArgumentError Weights(NaN, true, 0.0)   # positional form goes through the same checks
    @test_throws ArgumentError Weights(2.0, false, -1.0)
    # A vector of bases used to be accepted by the type and then hit a
    # MethodError inside ψ₂, ψ∞ and update!, whatever `updatew` was. It has no
    # defined meaning, so it is refused up front with a message that says so.
    @test_throws ArgumentError Weights(w=[2.0, 3.0])
    @test_throws ArgumentError Weights([2.0, 3.0], false, 1.0)
    @test Weights(w=2).w == 2                            # any Real scalar is fine
end

@testset "complex states through the stack" begin
    # The thesis's flagship runs (Kuramoto–Sivashinsky) push COMPLEX spectral
    # states through Parareal in anger; this pins that support in CI on a
    # complex Dahlquist. Explicit solvers both sides: no Jacobian machinery.
    λ = -0.4 + 2.0im
    zproblem = Dahlquist([1.0 + 0.0im], (0.0, 2.0); λ)
    zfine = solve(zproblem, finesolver)
    zparareal = Parareal(finesolver, coarsesolver;
                         parameters=PararealParameters(N=4, K=4), tolerance=Tolerance(ϵ=1e-12))
    zsolution = solve(zproblem, zparareal)
    u, t = flatten(zsolution)
    @test eltype(first(u)) == ComplexF64            # nothing silently realified
    @test norm(u[end] - zfine.u[end]) < 1e-8        # lands on the fine solve
    @test seams(zsolution) isa Vector{Float64}      # seams are norms: real, even for ℂ states
    @test all(≥(0), seams(zsolution))
    zthreads = solve(zproblem, zparareal; mode="THREADS")
    for n = 1:4
        @test zthreads[n].u == zsolution[n].u       # bitwise backend agreement holds for ℂ too
    end
end

@testset "update! on exact exponential orbit data (documented bias)" begin
    # On the EXACT orbit of u' = λu — U on the orbit, F its one-chunk fine
    # propagation — the estimator's ratio ‖F[i+1] − F[i]‖ / ‖U[i] − U[i−1]‖
    # spans TWO chunk gains (the F-difference sits one chunk downstream AND
    # one index up), so it returns the SQUARE of the one-chunk gain:
    # per-time w = e^{2λ}, not e^{λ}. Pinned here so any change to the
    # estimator is DELIBERATE, not accidental. Note the regime: smooth orbit
    # data. Mid-iteration Parareal data — iterate discrepancies, not orbit
    # spacings — is what the estimator actually sees in production, where a
    # perturbation-quotient reading applies instead.
    λ = 0.6
    T = [0.0, 1.0, 2.0, 3.0]
    U = [[exp(λ * T[n])] for n = 1:4]           # exact orbit boundary values
    F = [[exp(λ * (T[n] + 1.0))] for n = 1:4]   # exact one-chunk propagation of each
    w = Weights(updatew=true)
    NSDETimeParallel.update!(w, U, F, T)
    @test w.w ≈ exp(2λ) rtol=1e-10
end

@testset "superlinear contraction on a linear problem" begin
    # Classical Parareal theory on u' = λu: iterate errors obey a FACTORIAL
    # envelope e_k ≲ (Cᵏ/k!) e₀, not a geometric one — "converges eventually"
    # tests pass with a subtly wrong correction; a rate bound does not.
    # Self-calibrated: C is the first observed ratio, the slack is generous,
    # and the roundoff floor is excluded. A deliberately bad coarse solver
    # (one Euler step per chunk) keeps several iterations observable.
    zp = Dahlquist([1.0], (0.0, 3.0); λ=1.0)
    lp = Parareal(finesolver, Euler(h=0.5);
                  parameters=PararealParameters(N=6, K=6), tolerance=Tolerance(ϵ=0.0))
    e = solve(zp, lp).errors
    usable = [e[k] for k in eachindex(e) if isfinite(e[k]) && e[k] > 1e-11]
    @test length(usable) ≥ 3                        # enough iterations to see a rate
    C = usable[2] / usable[1]
    for k = 3:length(usable)
        @test usable[k] ≤ 100 * usable[1] * C^(k - 1) / factorial(k - 1) # factorial, not geometric
    end
    @test usable[end] / usable[end-1] < usable[2] / usable[1] # contraction improves with k
end

@testset "flatten" begin
    N = 4
    solution = solve(problem, make_parareal(N=N))
    u, t = flatten(solution)
    @test length(u) == length(t)
    @test length(t) == sum(length(solution[n].t) for n = 1:N) - (N - 1) # boundaries deduplicated
    @test issorted(t)
    @test u[1] == problem.u0
    @test norm(u[end] - finetruth(t[end])) < 1e-8
end

@testset "update! converts per-chunk ratios to per-time" begin
    w = Weights(updatew=true)
    U = [[0.0], [1.0], [2.0], [3.0]]
    F = [[0.0], [1.0], [3.0], [7.0]]  # per-chunk ratios: 2 (i=2) and 4 (i=3)
    T = [0.0, 2.0, 4.0, 6.0]          # chunks are 2 time units long
    NSDETimeParallel.update!(w, U, F, T)
    @test w.w ≈ 4^(1/2)               # max per-chunk ratio 4 → per-time base 2
    w1 = Weights(updatew=true)
    NSDETimeParallel.update!(w1, U, F, [0.0, 1.0, 2.0, 3.0])
    @test w1.w ≈ 4.0                  # unit chunks: conversion is the identity
end

@testset "cache and solution reuse" begin
    # make_parareal's default ϵ converges before K, so `errors` gets resized
    # down — a second run on the same solution used to throw BoundsError.
    parareal = make_parareal(N=5)
    fresh = solve(problem, parareal)
    cache = initialize_cache(problem, parareal)
    reused = initialize_solution(problem, parareal)
    parareal(cache, reused, problem)
    parareal(cache, reused, problem)   # the regression: reuse after down-resize
    @test reused.errors == fresh.errors # bitwise, run-to-run and vs fresh
    for n = 1:5
        @test reused[n].u == fresh[n].u
        @test reused[n].t == fresh[n].t
    end
    # saveiterates path: history capacity must also be restored on reuse
    withiter = initialize_solution(problem, parareal; saveiterates=true)
    parareal(cache, withiter, problem; saveiterates=true)
    k₁ = numiterates(withiter)
    parareal(cache, withiter, problem; saveiterates=true)
    @test numiterates(withiter) == k₁
    @test boundaryvalues(withiter.iterates[k₁]) == boundaryvalues(withiter.lastiterate)
end

@testset "convergence to the fine solve (serial)" begin
    N = 5
    parareal = make_parareal(N=N, K=N, ϵ=1e-10)
    solution = solve(problem, parareal)
    T = chunkgrid(problem, N)
    for t in T
        @test norm(solution(t) - finetruth(t)) < 1e-8
    end
    @test solution.errors[end] < solution.errors[1] # the iteration actually gains
    @test length(solution.errors) ≤ N
end

@testset "finite termination: k = N reproduces the fine solve" begin
    N = 4
    parareal = make_parareal(N=N, K=N, ϵ=1e-30, coarse=Euler(h=0.25)) # force full sweeps
    solution = solve(problem, parareal)
    for t in chunkgrid(problem, N)
        @test norm(solution(t) - finetruth(t)) ≤ 1e-10 # exact up to chunked compensated-sum residue
    end
end

@testset "chunk solves end on the chunk boundary bit for bit" begin
    # The exact-seam claim rests on this: the fine chunk endpoint must be the
    # final node, not an interpolation across a one-ulp overshoot.
    for n = 1:4
        T = chunkgrid(problem, 4)
        chunk = Logistic(0.3, (T[n], T[n+1]))
        sol = solve(chunk, finesolver)
        @test sol.t[end] == T[n+1]
        @test sol(T[n+1]) === sol.u[end]
    end
end

@testset "correction: exact prefix is assigned, never recomputed" begin
    # The interface freed by sweep k is copied from the fine output. Run with
    # ϵ = 0 and inspect the cache after each sweep through `saveiterates`:
    # every boundary value in the exact prefix must equal the fine endpoint
    # to the bit, sweep after sweep, with no ulp drift.
    N = 4
    parareal = make_parareal(N=N, K=N, ϵ=0.0)
    solution = parareal(problem; saveiterates=true)
    for k = 1:N
        iterate = solution.iterates[k]
        # after fine batch k, chunks 1..k are exact: their starts equal the previous chunk's fine endpoint
        for n = 2:k
            @test iterate[n].u[begin] == iterate[n-1].u[end]
        end
    end
    @test maxseam(solution) == 0.0
end

# Unwrap a NewtonFailure raised inside a threaded or distributed fine map,
# which the runtime may hand back wrapped once or twice.
function newtonfailure(err)
    err isa NewtonFailure && return err
    err isa TaskFailedException && return newtonfailure(err.task.result)
    err isa CompositeException && return newtonfailure(first(err))
    err isa RemoteException && return newtonfailure(err.captured.ex)
    err isa CapturedException && return newtonfailure(err.ex)
    return nothing
end

@testset "implicit solvers inside Parareal: success, and failure reaching the caller" begin
    # Success: backward Euler as the coarse propagator. Nothing in the
    # time-parallel layer assumes an explicit coarse solver.
    implicit = Parareal(finesolver, BackwardEuler(h=5e-2);
                        parameters=PararealParameters(N=5, K=5), tolerance=Tolerance(ϵ=1e-10))
    solution = solve(problem, implicit)
    for t in 0.0:0.1:1.0
        @test norm(solution(t) - finetruth(t)) < 1e-8
    end
    # Failure must surface, not be absorbed into a garbage chunk. u′ = u²
    # from u₀ = 1 has u(0.3) = 1/0.7; a backward-Euler step of 0.3 from u = 1
    # has no real stage (y = 1 + 0.3 y² needs 1 − 1.2 u ≥ 0), so the very
    # first coarse pass throws.
    blow = IVP((u, t) -> u .^ 2, [1.0], (0.0, 0.9))
    coarsefail = Parareal(RK4(h=1e-3), BackwardEuler(h=0.3);
                          parameters=PararealParameters(N=3, K=3), tolerance=Tolerance(ϵ=1e-10))
    @test_throws NewtonFailure solve(blow, coarsefail)
    # The same failure in the FINE map, under the serial and the threaded
    # backend: it reaches the caller either bare or wrapped by the runtime.
    finefail = Parareal(BackwardEuler(h=0.3), Euler(h=0.3);
                        parameters=PararealParameters(N=3, K=3), tolerance=Tolerance(ϵ=1e-10))
    @test_throws NewtonFailure solve(blow, finefail)
    err = try; solve(blow, finefail; mode="THREADS"); nothing; catch e; e; end
    @test newtonfailure(err) isa NewtonFailure
end

@testset "K cap without convergence behaves" begin
    parareal = make_parareal(N=6, K=2, ϵ=1e-30, coarse=Euler(h=0.25))
    solution = solve(problem, parareal)
    @test length(solution.errors) == 2
    @test all(isfinite, solution.errors)
    @test solution(1.0) isa AbstractVector # still a usable solution
end

@testset "backend agreement: THREADS ≡ SERIAL, bitwise" begin
    # One skeleton, disjoint chunk slots, serial correction: the float ops are
    # identical by construction, so this is equality, not approximation.
    N = 5
    sol_s = solve(problem, make_parareal(N=N); mode="SERIAL")
    sol_t = solve(problem, make_parareal(N=N); mode="THREADS")
    @test sol_s.errors == sol_t.errors
    for n = 1:N
        @test sol_s[n].u == sol_t[n].u
        @test sol_s[n].t == sol_t[n].t
    end
end

@testset "backend agreement: DISTRIBUTED ≡ SERIAL, worker-count independent" begin
    procs_added = addprocs(2; exeflags="--project=$(Base.active_project())")
    try
        @everywhere using NSDETimeParallel, NSDERungeKutta
        N = 4 # 4 chunks on 2 workers: the schedule must not affect the answer
        sol_s = solve(problem, make_parareal(N=N); mode="SERIAL")
        sol_d = solve(problem, make_parareal(N=N); mode="DISTRIBUTED")
        @test sol_s.errors == sol_d.errors
        for n = 1:N
            @test sol_s[n].u == sol_d[n].u
        end
    finally
        rmprocs(procs_added)
    end
end

@testset "makeGs guards injected starts" begin
    parareal = make_parareal(N=4)
    cache = NSDETimeParallel.PararealCache(problem, parareal)
    marker = fill(123.456, length(problem.u0))
    copyto!(cache.U[3], marker)
    cache.makeGs[3] = false
    NSDETimeParallel.coarseinit!(cache, parareal)
    @test cache.U[3] == marker          # injected guess kept
    @test cache.U[2] == cache.G[2]      # un-guarded chunk seeded from the coarse run
    @test cache.U[1] == problem.u0      # first chunk start copies u0 (no alias: MoWi writes into U[1])
    @test cache.U[1] !== problem.u0
end

@testset "saveiterates records the iteration history" begin
    parareal = make_parareal(N=5, K=5, ϵ=1e-12, coarse=Euler(h=0.2)) # several iterations
    solution = solve(problem, parareal; saveiterates=true)
    K = numiterates(solution)
    @test K ≥ 2
    @test length(solution.iterates) == K
    for n = 1:numchunks(solution)
        @test solution.iterates[K][n].u == solution.lastiterate[n].u # last snapshot = final answer
    end
    @test any(n -> solution.iterates[1][n].u != solution.iterates[K][n].u, 1:numchunks(solution))
end

@testset "collect_iterates! round-trips the opt-in disk history" begin
    parareal = make_parareal(N=3, K=3, ϵ=1e-30, coarse=Euler(h=0.25))
    reference = solve(problem, parareal; saveiterates=true)
    K = numiterates(reference)
    dir = mktempdir()
    for k = 1:K, n = 1:numchunks(reference)
        open(f -> serialize(f, (chunk_n = reference.iterates[k][n],)),
             joinpath(dir, "iter_$(k)_chunk_$(n).jls"), "w")
    end
    fresh = PararealSolution(problem, parareal; saveiterates=true)
    resize!(fresh.errors, K)
    resize!(fresh.iterates, K)
    collect_iterates!(fresh; directory=dir)
    for k = 1:K, n = 1:numchunks(reference)
        @test fresh.iterates[k][n].u == reference.iterates[k][n].u
    end
    @test_throws ArgumentError collect_iterates!(PararealSolution(problem, parareal); directory=dir)
end

@testset "cache reuse across a window shift (the MoWi seam)" begin
    N = 4
    parareal = make_parareal(N=N)
    injectstarts!(cache, source, τJ) = begin
        for n = 1:N
            tₙ = cache.T[n]
            if tₙ ≤ τJ
                copyto!(cache.U[n], source(tₙ))
                cache.makeGs[n] = false
            else
                cache.makeGs[n] = true
            end
        end
    end
    # window 1 on [0, 0.6], then shift the SAME cache to [0.3, 0.9]:
    p1 = copy(problem, problem.u0, 0.0, 0.6)
    cache = NSDETimeParallel.PararealCache(p1, parareal)
    w1 = parareal(cache, p1)
    p2 = copy(problem, problem.u0, 0.3, 0.9)
    NSDETimeParallel.shiftwindow!(cache, 0.3, 0.9)
    injectstarts!(cache, w1, 0.6)
    w2_reused = parareal(cache, p2)
    # reference: a FRESH cache with the same injected starts must agree bitwise
    fresh = NSDETimeParallel.PararealCache(p2, parareal)
    injectstarts!(fresh, w1, 0.6)
    w2_fresh = parareal(fresh, p2)
    @test w2_reused.errors == w2_fresh.errors
    for n = 1:N
        @test w2_reused[n].u == w2_fresh[n].u
        @test w2_reused[n].t == w2_fresh[n].t
    end
end

@testset "concrete types at the seams" begin
    parareal = make_parareal(N=4)
    solution = solve(problem, parareal)
    @test isconcretetype(eltype(solution.lastiterate.chunks)) # no Vector{Abstract…} boxes
    @test eltype(solution.errors) == Float64
    cache = NSDETimeParallel.PararealCache(problem, parareal)
    @test isconcretetype(eltype(cache.chunkproblems))
    @test isconcretetype(eltype(cache.finecaches))
end

@testset "boundary accessors" begin
    N = 5
    solution = solve(problem, make_parareal(N=N))
    T = boundarytimes(solution.lastiterate)
    U = boundaryvalues(solution.lastiterate)
    @test length(T) == length(U) == N + 1
    @test T ≈ chunkgrid(problem, N)                    # the Parareal grid, exactly
    @test all(u -> u isa AbstractVector, U)            # states, never spliced scalars
    @test U[1] == problem.u0                           # first boundary is u0
    for n in eachindex(T)
        @test norm(U[n] - finetruth(T[n])) < 1e-8      # converged boundaries sit on the fine solve
    end
    # Vector-valued problem: length must be N + 1, NOT N + d — vcat splices a
    # trailing state vector into its scalar components unless wrapped.
    vecproblem = SimplePendulum([π/4, 0.0], (0.0, 1.0))
    vecsolution = solve(vecproblem, Parareal(finesolver, coarsesolver;
        parameters=PararealParameters(N=3, K=3), tolerance=Tolerance(ϵ=1e-10)))
    V = boundaryvalues(vecsolution.lastiterate)
    @test length(V) == 4
    @test all(v -> v isa AbstractVector && length(v) == 2, V)
end

@testset "interpolation and edges" begin
    solution = solve(problem, make_parareal(N=5))
    @test solution(0.0) ≈ problem.u0
    @test norm(solution(1.0) - finetruth(1.0)) < 1e-8
    @test norm(solution(0.37) - finetruth(0.37)) < 1e-6
    @test solution(-5.0) == solution.lastiterate[1](-5.0)  # clamped below
    @test solution(99.0) ≈ solution(1.0)                   # clamped above
    @test Wnorm(solution.lastiterate, finetruth, 1.0) < 1e-6
end

@testset "loud failures" begin
    @test_throws ArgumentError solve(problem, make_parareal(N=3); mode="TYPO")
    # saveiterates is an argument error BEFORE any MPI check: it is invalid on
    # the pipelined backend regardless of environment.
    @test_throws ArgumentError solve(problem, make_parareal(N=3);
                                     backend=PipelinedMPIBackend(), saveiterates=true)
    if !NSDETimeParallel.MPI.Initialized()
        @test_throws ErrorException solve(problem, make_parareal(N=3); backend=MPIBackend())
        @test_throws ErrorException solve(problem, make_parareal(N=3); backend=PipelinedMPIBackend())
        @test_throws ErrorException solve(problem, make_parareal(N=3); mode="PIPELINED")
    end
end

@testset "ψ₁ guards zero states" begin
    Z = [zeros(2), ones(2)]
    fake = (U=Z, T=[0.0, 0.5, 1.0], U_=[zeros(2), zeros(2)], F=[zeros(2), ones(2)])
    @test isfinite(ψ₁(fake, 1, Weights()))
end

@testset "ψ∞: sandwich, weighting, convergence" begin
    # Exact-arithmetic sandwich on a synthetic cache: ψ₂ ≤ ψ∞ ≤ N·ψ₂, and the
    # max picks the boundary the mean averages away.
    U = [[1.0], [2.0], [3.0]]
    F = [[1.0], [2.5], [3.0]]                       # one defect of 0.5 at n = 2
    fake = (U=U, T=[0.0, 1.0, 2.0, 3.0], F=F)
    unweighted = Weights()                          # w = 1
    @test ψ∞(fake, 1, unweighted) ≈ 0.5             # the max IS the lone defect
    @test ψ₂(fake, 1, unweighted) ≈ 0.5 / 3         # the mean dilutes it by N
    @test ψ₂(fake, 1, unweighted) ≤ ψ∞(fake, 1, unweighted) ≤ 3 * ψ₂(fake, 1, unweighted)
    weighted = Weights(w=exp(1.0))
    @test ψ∞(fake, 1, weighted) ≈ exp(-1.0) * 0.5   # weight w^(T[1] − T[2]) at n = 2
    @test ψ∞(fake, 1, Weights(w=0.5)) ≈ ψ∞(fake, 1, unweighted) # sub-unit base clamped to 1

    # On a real run, ψ∞ converges to the fine solve and, at the FIRST sweep
    # (both criteria measure the same iterate-0 defects), the sandwich holds
    # between the recorded errors of two otherwise identical runs.
    N = 5
    mk(ψ) = Parareal(finesolver, coarsesolver;
                     parameters=PararealParameters(N=N, K=N),
                     tolerance=Tolerance(ϵ=1e-10, ψ=ψ, weights=Weights(w=1.0)))
    s∞ = solve(problem, mk(ψ∞))
    s₂ = solve(problem, mk(ψ₂))
    @test s₂.errors[1] ≤ s∞.errors[1] ≤ N * s₂.errors[1] + eps()
    @test numiterates(s∞) ≥ numiterates(s₂)         # the max is the stricter rule at equal ϵ
    u, t = flatten(s∞)
    @test norm(u[end] - finetruth(t[end])) < 1e-8
end

@testset "RecipesBase recipes (headless)" begin
    psol = solve(problem, make_parareal(N=4))
    for obj in (psol, psol.lastiterate,
                NSDEBase._PhasePlot(psol), NSDEBase._PhasePlot(psol.lastiterate),
                NSDEBase._Convergence(psol))
        @test RecipesBase.apply_recipe(Dict{Symbol,Any}(), obj) isa Vector{RecipesBase.RecipeData}
    end
end


end # outer testset
