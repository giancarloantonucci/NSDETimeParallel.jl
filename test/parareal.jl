@testset "PararealParameters validation" begin
    @test_throws ArgumentError PararealParameters(N=0)          # default K = N = 0
    @test_throws ArgumentError PararealParameters(N=5, K=0)
    @test_throws ArgumentError PararealParameters(-1, 3)
    @test_logs (:warn, r"dead work") PararealParameters(N=2, K=5)
    @test PararealParameters(N=5, K=5) isa PararealParameters   # K = N: silent
end

@testset "complex states through the stack" begin
    # Kuramoto–Sivashinsky discretisations push COMPLEX spectral
    # states through Parareal in anger; this pins that support in CI on a
    # complex Dahlquist. Explicit solvers both sides: no Jacobian machinery.
    λ = -0.4 + 2.0im
    zproblem = Dahlquist([1.0 + 0.0im], (0.0, 2.0); λ)
    zfine = solve(zproblem, finesolver)
    zparareal = Parareal(
        finesolver,
        coarsesolver;
        parameters=PararealParameters(N=4, K=4),
        criterion=MeanRelativeUpdate(ϵ=1e-12),
    )
    zsolution = solve(zproblem, zparareal)
    u, t = flatten(zsolution)
    @test eltype(first(u)) == ComplexF64            # nothing silently realified
    @test norm(u[end] - zfine.u[end]) < 1e-8        # lands on the fine solve
    @test seams(zsolution) isa Vector{Float64}      # seams are norms: real, even for ℂ states
    @test all(≥(0), seams(zsolution))
    zthreads = solve(zproblem, zparareal; backend=ThreadsBackend())
    for n in 1:4
        @test zthreads[n].u == zsolution[n].u       # bitwise backend agreement holds for ℂ too
    end
end

@testset "superlinear contraction on a linear problem" begin
    # Classical Parareal theory on u' = λu: iterate errors obey a FACTORIAL
    # envelope e_k ≲ (Cᵏ/k!) e₀, not a geometric one — "converges eventually"
    # tests pass with a subtly wrong correction; a rate bound does not.
    # Self-calibrated: C is the first observed ratio, the slack is generous,
    # and the roundoff floor is excluded. A deliberately bad coarse solver
    # (one Euler step per chunk) keeps several iterations observable.
    zp = Dahlquist([1.0], (0.0, 3.0); λ=1.0)
    lp = Parareal(
        finesolver,
        Euler(h=0.5);
        parameters=PararealParameters(N=6, K=6),
        criterion=MeanRelativeUpdate(ϵ=0.0),
    )
    e = solve(zp, lp).errors
    usable = [e[k] for k in eachindex(e) if isfinite(e[k]) && e[k] > 1e-11]
    @test length(usable) ≥ 3                        # enough iterations to see a rate
    C = usable[2] / usable[1]
    for k in 3:length(usable)
        @test usable[k] ≤ 100 * usable[1] * C^(k - 1) / factorial(k - 1) # factorial, not geometric
    end
    @test usable[end] / usable[end - 1] < usable[2] / usable[1] # contraction improves with k
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

@testset "correction: exact prefix is assigned, never recomputed" begin
    # The interface freed by sweep k is copied from the fine output. Run with
    # ϵ = 0 and inspect the cache after each sweep through `saveiterates`:
    # every boundary value in the exact prefix must equal the fine endpoint
    # to the bit, sweep after sweep, with no ulp drift.
    N = 4
    parareal = make_parareal(N=N, K=N, criterion=NoStop())
    solution = parareal(problem; history=:iterates)
    for k in 1:N
        iterate = solution.iterates[k]
        # after fine batch k, chunks 1..k are exact: their starts equal the previous chunk's fine endpoint
        for n in 2:k
            @test iterate[n].u[begin] == iterate[n - 1].u[end]
        end
    end
    @test maxseam(solution) == 0.0
end

# Unwrap a NewtonFailure raised inside a threaded or distributed fine map,
# which the runtime may hand back wrapped once or twice.

@testset "implicit solvers inside Parareal: success and explicit failure status" begin
    # Success: backward Euler as the coarse propagator. Nothing in the
    # time-parallel layer assumes an explicit coarse solver.
    implicit = Parareal(
        finesolver,
        BackwardEuler(h=5e-2);
        parameters=PararealParameters(N=5, K=5),
        criterion=MeanRelativeUpdate(ϵ=1e-10),
    )
    solution = solve(problem, implicit)
    for t in 0.0:0.1:1.0
        @test norm(solution(t) - finetruth(t)) < 1e-8
    end
    # Failure must surface, not be absorbed into a garbage chunk. u′ = u²
    # from u₀ = 1 has u(0.3) = 1/0.7; a backward-Euler step of 0.3 from u = 1
    # has no real stage (y = 1 + 0.3 y² needs 1 − 1.2 u ≥ 0), so the very
    # first coarse pass fails.
    blow = IVP((u, t) -> u .^ 2, [1.0], (0.0, 0.9))
    coarsefail = Parareal(
        RK4(h=1e-3),
        BackwardEuler(h=0.3);
        parameters=PararealParameters(N=3, K=3),
        criterion=MeanRelativeUpdate(ϵ=1e-10),
    )
    failed = @test_logs (:warn, r"Chunk solve failed") match_mode=:any solve(blow, coarsefail)
    @test failed.diagnostics.status == :nonfinite
    # The same failure in the FINE map, under the serial and the threaded
    # backend: it returns the same explicit nonfinite status.
    finefail = Parareal(
        BackwardEuler(h=0.3),
        Euler(h=0.3);
        parameters=PararealParameters(N=3, K=3),
        criterion=MeanRelativeUpdate(ϵ=1e-10),
    )
    failed = @test_logs (:warn, r"Chunk solve failed") match_mode=:any solve(blow, finefail)
    @test failed.diagnostics.status == :nonfinite
    failed = @test_logs (:warn, r"Chunk solve failed") match_mode=:any solve(
        blow, finefail; backend=ThreadsBackend()
    )
    @test failed.diagnostics.status == :nonfinite
end

@testset "K cap without convergence behaves" begin
    parareal = make_parareal(N=6, K=2, ϵ=1e-30, coarse=Euler(h=0.25))
    solution = solve(problem, parareal)
    @test length(solution.errors) == 2
    @test solution.errors[1] == Inf # no update exists before the first correction
    @test all(isfinite, solution.errors[2:end])
    @test solution.diagnostics.status == :budget_exhausted
    @test solution(1.0) isa AbstractVector # still a usable solution
end

@testset "loud failures" begin
    @test_throws ArgumentError solve(problem, make_parareal(N=3); history=:bad)
    # saveiterates is an argument error BEFORE any MPI check: it is invalid on
    # the pipelined backend regardless of environment.
    @test_throws ArgumentError solve(
        problem, make_parareal(N=3); backend=PipelinedMPIBackend(), history=:iterates
    )
    if !NSDETimeParallel.MPI.Initialized()
        @test_throws ErrorException solve(problem, make_parareal(N=3); backend=MPIBackend())
        @test_throws ErrorException solve(
            problem, make_parareal(N=3); backend=PipelinedMPIBackend()
        )
        @test_throws ErrorException solve(
            problem, make_parareal(N=3); backend=PipelinedMPIBackend()
        )
    end
end
