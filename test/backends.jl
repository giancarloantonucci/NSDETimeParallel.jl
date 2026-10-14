@testset "backend agreement: THREADS ≡ SERIAL, bitwise" begin
    # One skeleton, disjoint chunk slots, serial correction: the float ops are
    # identical by construction, so this is equality, not approximation.
    N = 5
    sol_s = solve(problem, make_parareal(N=N); backend=SerialBackend())
    sol_t = solve(problem, make_parareal(N=N); backend=ThreadsBackend())
    @test sol_s.errors == sol_t.errors
    for n in 1:N
        @test sol_s[n].u == sol_t[n].u
        @test sol_s[n].t == sol_t[n].t
    end
end

@testset "backend agreement: DISTRIBUTED ≡ SERIAL, worker-count independent" begin
    procs_added = addprocs(2; exeflags="--project=$(Base.active_project())")
    try
        @everywhere using NSDETimeParallel, NSDERungeKutta
        N = 4 # 4 chunks on 2 workers: the schedule must not affect the answer
        sol_s = solve(problem, make_parareal(N=N); backend=SerialBackend())
        sol_d = solve(problem, make_parareal(N=N); backend=DistributedBackend())
        @test sol_s.errors == sol_d.errors
        base = IVP(-1.0, [1.0], (0, 2))
        solver = Parareal(RK4(h=0.03), Euler(h=0.1); N, criterion=NoStop())
        cache = NSDEBase.initialize_cache(base, solver)
        changed = IVP(-2.0, [0.4], (0, 2))
        @test boundaryvalues(solver(cache, changed; backend=DistributedBackend()).lastiterate) ==
            reference(changed, solver)
        for n in 1:N
            @test sol_s[n].u == sol_d[n].u
        end
        pbad=IVP((u, t)->u .^ 2, [1.0], (0.0, 0.4))
        badsolver=Parareal(
            BackwardEuler(h=0.1, Mₙ=1, εᵣ=1e-14), Euler(h=0.05); N=4, criterion=NoStop()
        )
        failed=solve(pbad, badsolver; backend=DistributedBackend())
        bug=IVP((u, t)->throw(BoundsError(u, 99)), [1.0], (0.0, 0.1))
        probe=Parareal(RK4(h=0.01), Euler(h=0.02); N=1, criterion=NoStop())
        err = try
            solve(bug, probe; backend=DistributedBackend())
            nothing
        catch e
            e
        end
        @test err isa RemoteException
        @test err.captured.ex isa BoundsError
        @test failed.diagnostics.status == :nonfinite
        @test all(u->all(isnan, u), boundaryvalues(failed.lastiterate))
    finally
        rmprocs(procs_added)
    end
end

@testset "Adaptive state is local to a solve and thread" begin
    fine = DormandPrince54(h=0.1, εᵣ=1e-7, save_stepsizes=true)
    p = Logistic(0.2, (0.0, 2.0))
    a, b = solve(p, fine), solve(p, fine)
    @test a.t == b.t
    @test a.u == b.u
    @test fine.stepsize.h == 0.1
    @test a.stepsizes !== b.stepsizes
    @test a.t[end] == 2.0
    s = Parareal(fine, RK4(h=0.1); N=8, criterion=NoStop())
    x = solve(p, s)
    y = solve(p, s; backend=ThreadsBackend())
    @test boundaryvalues(x.lastiterate) == boundaryvalues(y.lastiterate) == reference(p, s)
    @test fine.stepsize.h == 0.1
end

@testset "CI exercises multiple Julia threads" begin
    if get(ENV, "NSDE_REQUIRE_THREADS", "false")=="true"
        @test Threads.nthreads() >= 2
    end
end
