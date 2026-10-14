@testset "chunk solves end on the chunk boundary bit for bit" begin
    # The exact-seam claim rests on this: the fine chunk endpoint must be the
    # final node, not an interpolation across a one-ulp overshoot.
    for n in 1:4
        T = chunkgrid(problem, 4)
        chunk = Logistic(0.3, (T[n], T[n + 1]))
        sol = solve(chunk, finesolver)
        @test sol.t[end] == T[n + 1]
        @test sol(T[n + 1]) === sol.u[end]
    end
end

@testset "Automatic chunk plans and exact reference" begin
    for dt in (0.23, 1.0), N in (1, 4), start in (0.0, 103.0)
        p = Logistic([0.2], (start, start + N*dt))
        s = Parareal(RK4(h=0.03), RK4(h=0.1); N, criterion=NoStop())
        a = solve(p, s; history=:interfaces)
        b = solve(p, s; backend=ThreadsBackend())
        @test boundaryvalues(a.lastiterate) == reference(p, s)
        @test boundaryvalues(b.lastiterate) == boundaryvalues(a.lastiterate)
        @test a.diagnostics.status == :exact
        @test a.diagnostics.fine_batches == N
        @test a.diagnostics.sweeps == N-1
        @test a.diagnostics.coarse_chunk_calls == N*(N-1)÷2
        @test a.errors === a.diagnostics.residuals
        @test length(a.diagnostics.interfaces) == N
        @test a.lastiterate[end].t[end] == p.tspan[end]
    end
end

@testset "Reference preserves absolute time" begin
    rhs = NSDEBase.RHS((u, t)->[sin(t)+u[1]])
    p = NSDEBase.IVP(rhs, [0.2], (3.0, 4.0))
    s = Parareal(RK4(h=0.03), RK4(h=0.1); N=4, criterion=NoStop())
    @test boundaryvalues(s(p).lastiterate) == reference(p, s)
    @test norm(reference(p, s)[end] - reference(NSDEBase.remake(p; tspan=(0.0, 1.0)), s)[end]) > 0.1
end

@testset "Scheduled chunks in native Parareal" begin
    for dt in (1.6, 500/99), N in (3, 7), coarse in (RK4(h=0.1), Euler(h=0.01))
        F=chunk_solver(RK4(h=0.001), dt)
        G=chunk_solver(coarse, dt)
        p=Lorenz([2.0, 3.0, -14.0], (0.0, N*dt))
        a=Parareal(F, G; parameters=PararealParameters(N=N), criterion=NoStop())
        s=solve(p, a; backend=ThreadsBackend(), history=:interfaces)
        u=[2.0, 3.0, -14.0]
        for n in 1:N
            u=solve(Lorenz(u, (0.0, dt)), F).u[end]
        end
        @test boundaryvalues(s.lastiterate)[end]==u
        @test s.diagnostics.status==:exact
        @test maxseam(s.lastiterate)==0
    end
    for method in (BackwardEuler(h=0.01), RadauIIA3(h=0.01))
        c=chunk_solver(method, 0.123)
        a=solve(Logistic([0.1], (0.0, 0.123)), c)
        b=solve(Logistic([0.1], (10.0, 10.123)), c)
        @test a.u[end]==b.u[end]
        @test a.t[end]==0.123
    end
end

@testset "Nominal plans survive 20,000 shifted windows" begin
    p=Logistic(0.3, (0.0, 0.6))
    solver=Parareal(RK4(h=1e-4), RK4(h=1e-3); N=10, criterion=NoStop())
    cache=initialize_cache(p, solver)
    finecache=cache.finecaches[1]
    oldplan=cache.finesolver
    unchanged=true
    spans_valid=true
    for i in 0:19999
        t0=i*0.3
        shiftwindow!(cache, t0, t0+0.6; duration=0.6/10)
        NSDETimeParallel.prepare_chunks!(cache, solver)
        unchanged &= cache.finecaches[1] === finecache && cache.finesolver === oldplan
        spans_valid &= cache.finesolver.plan.full == 600 && cache.finesolver.plan.tail == 0
        if i % 199 == 0
            chunk=cache.chunkproblems[1]
            sol=solve(chunk, cache.finesolver)
            spans_valid &= length(sol.t)==601 && sol.t[end]==chunk.tspan[2]
        end
    end
    @test unchanged
    @test spans_valid
    @test cache.duration == 0.6/10
    shiftwindow!(cache, 6000.0, 6000.3; duration=0.3/10)
    NSDETimeParallel.prepare_chunks!(cache, solver)
    @test cache.finecaches[1] !== finecache
    @test cache.finesolver.plan.full == 300
    @test_throws ArgumentError shiftwindow!(cache, 6000.0, 6000.6)
end
