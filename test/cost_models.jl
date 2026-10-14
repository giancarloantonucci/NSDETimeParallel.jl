@testset "count-based speed-up models" begin
    @test theoretical_speedup(1, 10, 0.0) == 10.0
    # Two fine batches: 9 initial coarse calls, then 8 correction calls.
    @test theoretical_speedup(2, 10, 0.1) ≈ 10 / (2 + 17*0.1)
    @test theoretical_speedup(2, 10, 0.1; estimate=:work) ≈ 10 / (19/10 + 17*0.1)
    @test theoretical_speedup(1, 1, 7.0) == 1.0 # no coarse call
    @test theoretical_speedup(10, 10, 0.0) == 1.0 # exact completion
    @test theoretical_speedup(10, 10, 0.1) ≈ 10/(10+45*0.1)
    @test theoretical_speedup(2, 10, 0.1; wallclock=false) ==
        theoretical_speedup(2, 10, 0.1; estimate=:work)
    @test theoretical_speedup(2, 10, 0.1; wallclock=false, estimate=:wallclock) ==
        theoretical_speedup(2, 10, 0.1)
    for args in
        ((0, 10, 0.1), (11, 10, 0.1), (1, 0, 0.1), (1, 10, -1.0), (1, 10, NaN), (1, 10, Inf))
        @test_throws ArgumentError theoretical_speedup(args...)
    end
    @test_throws ArgumentError theoretical_speedup(2, 10, 0.1; estimate=:aubanel)
end

@testset "model agrees with actual recorded calls" begin
    for N in (1, 2, 5), K in 1:N
        p = Logistic(0.3, (0.0, 1.0))
        solver = Parareal(RK4(h=0.01), Euler(h=0.05); N, K, criterion=NoStop())
        sol = solve(p, solver)
        d = sol.diagnostics
        @test d.fine_batches == K
        @test d.sweeps == K-1
        # Enumerate the initial pass and the actual correction loop ranges.
        coarse_calls = (N-1) + sum((length((k + 1):(N - 1)) for k in 1:(K - 1)); init=0)
        fine_calls = sum(length(k:N) for k in 1:K)
        @test d.coarse_chunk_calls == coarse_calls
        @test d.fine_chunk_calls == fine_calls
        for ζ in (0.0, 0.1, 1.0), model in (:wallclock, :work)
            expected = N/((model === :wallclock ? K : fine_calls/N) + ζ*coarse_calls)
            @test theoretical_speedup(K, N, ζ; estimate=model) == expected
            @test theoretical_speedup(d, ζ; estimate=model) == expected
            @test theoretical_speedup(sol, ζ; estimate=model) == expected
        end
    end
    # A weighted rule may accept the coarse starts in its first fine batch.
    solver = Parareal(RK4(h=0.01), Euler(h=0.05); N=4, criterion=WeightedDefect(ϵ=1.0))
    sol = solve(problem, solver)
    @test sol.diagnostics.sweeps == 0
    @test theoretical_speedup(sol, 0.1) == 4/(1+3*0.1)
    unstarted = initialize_solution(problem, solver)
    @test_throws ArgumentError theoretical_speedup(unstarted, 0.1)
    failed = deepcopy(sol.diagnostics)
    failed.status = :nonfinite
    @test_throws ArgumentError theoretical_speedup(failed, 0.1)
    pipeline = deepcopy(sol.diagnostics)
    pipeline.fine_batches = nothing
    @test_throws ArgumentError theoretical_speedup(pipeline, 0.1)
    # The diagnostics overload uses recorded counts, not a recomputed formula.
    changed = deepcopy(sol.diagnostics)
    changed.coarse_chunk_calls += 1
    @test theoretical_speedup(changed, 0.1) == 4/(1+4*0.1)
end

@testset "costratio" begin
    ζ = costratio(problem, make_parareal(N=5))
    @test isfinite(ζ) && ζ > 0
    @test ζ < 1                        # the coarse solver takes 50× fewer steps here
    @test_throws ArgumentError costratio(problem, make_parareal(N=5); repeats=0)
end

@testset "fine-batch cost budget" begin
    # Exhaustive inversion against independently enumerated call counts.
    for N in (1, 2, 7, 64),
        ζ in (0.0, 0.02, 0.3), S in (0.05, 0.9, 1.0, 1.7, 4.0, 60.0),
        model in (:wallclock, :work)

        feasible = Int[]
        for K in 1:N
            coarse = N-1 + sum((N-k-1 for k in 1:(K - 1)); init=0)
            fine = sum(N-k+1 for k in 1:K)
            speedup = N/((model === :wallclock ? K : fine/N)+ζ*coarse)
            speedup >= S && push!(feasible, K)
        end
        expected = isempty(feasible) ? 0 : maximum(feasible)
        @test iteration_budget(S, N, ζ; estimate=model) == expected
    end
    @test iteration_budget(1e-300, 10, 0.1) == 10
    @test iteration_budget(1e6, 10, 0.1) == 0
    @test iteration_budget(2, 10, 0.1) == 2
    @test iteration_budget(1.0, 1, 1.0) == 1
    @test iteration_budget(2, 10, 0.1; estimate=:naive) == iteration_budget(2, 10, 0.1)
    @test iteration_budget(2, 10, 0.1; estimate=:improved) ==
        iteration_budget(2, 10, 0.1; estimate=:work)
    for S in (0.0, Inf, NaN)
        @test_throws ArgumentError iteration_budget(S, 10, 0.1)
    end
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
