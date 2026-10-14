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
    for n in 1:5
        @test reused[n].u == fresh[n].u
        @test reused[n].t == fresh[n].t
    end
    # saveiterates path: history capacity must also be restored on reuse
    withiter = initialize_solution(problem, parareal; history=:iterates)
    parareal(cache, withiter, problem; history=:iterates)
    k₁ = numiterates(withiter)
    parareal(cache, withiter, problem; history=:iterates)
    @test numiterates(withiter) == k₁
    @test boundaryvalues(withiter.iterates[k₁]) == boundaryvalues(withiter.lastiterate)
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
    solution = solve(problem, parareal; history=:iterates)
    K = numiterates(solution)
    @test K ≥ 2
    @test length(solution.iterates) == K
    for n in 1:numchunks(solution)
        @test solution.iterates[K][n].u == solution.lastiterate[n].u # last snapshot = final answer
    end
    @test any(n -> solution.iterates[1][n].u != solution.iterates[K][n].u, 1:numchunks(solution))
end

@testset "collect_iterates! round-trips the opt-in disk history" begin
    parareal = make_parareal(N=3, K=3, ϵ=1e-30, coarse=Euler(h=0.25))
    reference = solve(problem, parareal; history=:iterates)
    K = numiterates(reference)
    dir = mktempdir()
    for k in 1:K, n in 1:numchunks(reference)
        open(
            f -> serialize(f, (chunk_n=reference.iterates[k][n],)),
            joinpath(dir, "iter_$(k)_chunk_$(n).jls"),
            "w",
        )
    end
    fresh = PararealSolution(problem, parareal; history=:iterates)
    resize!(fresh.errors, K)
    resize!(fresh.iterates, K)
    collect_iterates!(fresh; directory=dir)
    for k in 1:K, n in 1:numchunks(reference)
        @test fresh.iterates[k][n].u == reference.iterates[k][n].u
    end
    @test_throws ArgumentError collect_iterates!(PararealSolution(problem, parareal); directory=dir)
end

@testset "cache reuse across a window shift (the MoWi seam)" begin
    N = 4
    parareal = make_parareal(N=N)
    injectstarts!(cache, source, τJ) = begin
        for n in 1:N
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
    p1 = remake(problem; u0=problem.u0, tspan=(0.0, 0.6))
    cache = NSDETimeParallel.PararealCache(p1, parareal)
    w1 = parareal(cache, p1)
    p2 = remake(problem; u0=problem.u0, tspan=(0.3, 0.9))
    NSDETimeParallel.shiftwindow!(cache, 0.3, 0.9)
    injectstarts!(cache, w1, 0.6)
    w2_reused = parareal(cache, p2)
    # reference: a FRESH cache with the same injected starts must agree bitwise
    fresh = NSDETimeParallel.PararealCache(p2, parareal)
    injectstarts!(fresh, w1, 0.6)
    w2_fresh = parareal(fresh, p2)
    @test w2_reused.errors == w2_fresh.errors
    for n in 1:N
        @test w2_reused[n].u == w2_fresh[n].u
        @test w2_reused[n].t == w2_fresh[n].t
    end
end

@testset "Weight state resets between runs" begin
    p = Lorenz([2.0, 3.0, -14.0], (0.0, 1.0))
    for weight in (SecantWeight(), NeighbourWeight(), FixedWeight(2.0))
        s = Parareal(
            RK4(h=0.002), RK4(h=0.02); N=8, criterion=WeightedDefect(ϵ=1e-10, weight=weight)
        )
        a = solve(p, s; history=:interfaces)
        b = solve(p, s; backend=ThreadsBackend(), history=:interfaces)
        @test a.diagnostics.weights == b.diagnostics.weights
        @test a.diagnostics.interfaces == b.diagnostics.interfaces
        @test s.criterion.weight.w == weight.w
        cache = NSDEBase.initialize_cache(p, s)
        c = s(cache, p; history=:interfaces)
        d = s(cache, p; history=:interfaces)
        @test c.diagnostics.weights == d.diagnostics.weights
        @test c.diagnostics.interfaces == d.diagnostics.interfaces
    end
end

@testset "cache reuse across criterion and solver types" begin
    p = Lorenz([2.0, 3.0, -14.0], (0.0, 0.4))
    probe = Parareal(
        RK4(h=0.002), Euler(h=0.02); N=4, criterion=WeightedDefect(ϵ=1e-12, weight=SecantWeight())
    )
    cache = initialize_cache(p, probe)
    a = probe(cache, p)
    frozen = remake(
        probe; criterion=WeightedDefect(ϵ=1e-12, weight=FixedWeight(a.diagnostics.weights[end]))
    )
    adaptive = remake(
        frozen;
        finesolver=DormandPrince54(h=0.001, εₐ=1e-12, εᵣ=1e-10),
        coarsesolver=DormandPrince54(h=0.01, εₐ=1e-10, εᵣ=1e-7),
    )
    buffer = initialize_solution(p, probe; history=:iterates)
    for solver in (frozen, adaptive, probe, remake(probe; criterion=RelativeUpdate()))
        fresh = solve(p, solver; history=:iterates)
        reused = solver(cache, buffer, p; history=:iterates)
        @test reused.errors == fresh.errors
        @test reused.diagnostics.weights == fresh.diagnostics.weights
        @test reused.diagnostics.status == fresh.diagnostics.status
        @test boundaryvalues(reused.lastiterate) == boundaryvalues(fresh.lastiterate)
        @test reused.iterates[end].chunks[end].u == fresh.iterates[end].chunks[end].u
    end
    @test_throws DimensionMismatch remake(probe; N=5)(cache, p)
end
