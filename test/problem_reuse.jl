@testset "Reused cache honours the supplied problem" begin
    p = IVP(-1.0, [1.0], (0, 2))
    solver = Parareal(RK4(h=0.03), Euler(h=0.1); N=4, criterion=NoStop())
    @test boundaryvalues(solve(p, solver).lastiterate) == reference(p, solver)
    @test isfinite(costratio(p, solver; repeats=1))
    cache = NSDEBase.initialize_cache(p, solver)
    result = NSDEBase.initialize_solution(p, solver)
    for backend in (SerialBackend(), ThreadsBackend(), DistributedBackend())
        for q in (remake(p; u0=[0.7]), IVP(-2.0, [0.4], (0, 2)))
            solver(cache, result, q; backend)
            @test boundaryvalues(result.lastiterate) == reference(q, solver)
            @test result[1].u[1] == q.u0
            @test q.u0 == (q.rhs === p.rhs ? [0.7] : [0.4])
        end
    end
    @test_throws ArgumentError solver(cache, result, remake(p; tspan=(1, 3)))
    @test_throws DimensionMismatch solver(cache, result, remake(p; u0=[1.0, 2.0]))
    @test_throws ArgumentError NSDETimeParallel.PararealCache(p, solver; duration=0.25)
    @test_throws ArgumentError NSDETimeParallel.PararealCache(p, solver; duration=Inf)
    @test remake(solver; N=7).parameters.K == 7
    @test remake(solver; N=7, K=2).parameters.K == 2
    @test remake(remake(solver; K=2); N=7).parameters.K == 2
    @test theoretical_speedup(result, 0.1; wallclock=false) ==
        theoretical_speedup(result, 0.1; estimate=:work)
    @test theoretical_speedup(result.diagnostics, 0.1; wallclock=false, estimate=:wallclock) ==
        theoretical_speedup(result, 0.1)
end

@testset "Nonfinite fine batch keeps residual and weight histories aligned" begin
    p = IVP((u, t)->fill(NaN, length(u)), [1.0], (0.0, 1.0))
    # One chunk bypasses the coarse pass: failure occurs in the first fine batch.
    solver = Parareal(RK4(h=0.1), Euler(h=0.1); N=1, criterion=WeightedDefect())
    result = solve(p, solver)
    @test result.diagnostics.status == :nonfinite
    @test result.errors == [Inf]
    @test length(result.diagnostics.weights) ==
        length(result.errors) ==
        result.diagnostics.fine_batches
    @test result.diagnostics.weights == [1.0]
end

@testset "RHS-dependent exponential caches are rebuilt" begin
    p = IVP(-1.0, [1.0], (0, 2))
    q = IVP(-2.0, [0.5], (0, 2))
    solver = Parareal(LawsonEuler(h=0.07), Euler(h=0.1); N=4, criterion=NoStop())
    cache = NSDEBase.initialize_cache(p, solver)
    solver(cache, p)
    result = solver(cache, q)
    @test boundaryvalues(result.lastiterate) == reference(q, solver)
    @test result[end].u[end][1] ≈ 0.5exp(-4)
end
