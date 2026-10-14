@testset "Criteria and backend contracts" begin
    p = Logistic(0.2, (0.0, 1.0))
    s = Parareal(RK4(h=0.01), RK4(h=0.1); N=4)
    @test s.criterion isa RelativeUpdate
    r=solve(p, s)
    @test r.diagnostics.sweeps >= 1
    @test r.diagnostics.fine_batches == r.diagnostics.sweeps+1
    @test_throws ArgumentError RelativeUpdate(atol=0, rtol=0)
    @test_throws ArgumentError WeightedDefect(p=0.5)
    @test_throws ArgumentError NSDETimeParallel.pipeline_kind(RelativeUpdate())
    @test_throws ArgumentError NSDETimeParallel.pipeline_kind(
        WeightedDefect(weight=SecantWeight(), p=Inf)
    )
    @test NSDETimeParallel.pipeline_kind(WeightedDefect(p=Inf)) == :defect
    @test_throws ArgumentError solve(p, s; history=:bad)
    v=solve(p, s; history=:iterates)
    @test length(v.iterates)==v.diagnostics.fine_batches
    @test Wnorm(v.lastiterate, reference(p, s), FixedWeight(2.0)) ==
        Wnorm(v.lastiterate, reference(p, s), 2.0)
end

@testset "Residual formulae and unchanged snapshots" begin
    U=[[1.0], [2.0], [3.0]]
    F=[[1.0], [2.5], [3.0]]
    cache=(
        U=U,
        F=F,
        U_=deepcopy(U),
        T=[0.0, 1.0, 2.0, 3.0],
        weightstate=NSDETimeParallel.WeightState(FixedWeight()),
    )
    avg=WeightedDefect()
    mx=WeightedDefect(p=Inf)
    @test residual(mx, cache, 1) ≈ 0.5
    @test residual(avg, cache, 1) ≈ 0.5/3
    @test residual(avg, cache, 1) ≤ residual(mx, cache, 1) ≤ 3residual(avg, cache, 1)
    cache.weightstate.w=exp(1.0)
    @test residual(mx, cache, 1) ≈ exp(-1.0)*0.5
    cache.weightstate.w=0.5
    @test residual(mx, cache, 1) ≈ 0.5
    @test residual(WeightedDefect(p=2), cache, 1) ≈ 0.5/sqrt(3)
    old=deepcopy(cache.U_)
    @test residual(MeanRelativeUpdate(), cache, 2) == 0
    @test cache.U_ == old
    fake=(U=[zeros(2), ones(2)], F=[zeros(2), ones(2)], U_=[zeros(2), zeros(2)])
    @test isfinite(residual(MeanRelativeUpdate(), fake, 1))
    fake=(U=[[0.0], [1.0], [1.0]], U_=[[0.0], [1.0], [0.7]])
    check=RelativeUpdate(atol=0, rtol=0.1)
    @test residual(check, fake, 1) == Inf
    @test residual(check, fake, 2) ≈ 3
    @test !accept(check, residual(check, fake, 2))
    @test residual(RelativeUpdate(atol=0, rtol=1), (U=[[0.0]], U_=[[0.0]]), 2) == 0
    mk(p) = Parareal(finesolver, coarsesolver; N=5, criterion=WeightedDefect(ϵ=1e-10, p=p))
    a=solve(problem, mk(1))
    b=solve(problem, mk(Inf))
    @test a.errors[1] ≤ b.errors[1] ≤ 5a.errors[1]+eps()
    @test numiterates(b) ≥ numiterates(a)
    u, t=flatten(b)
    @test norm(u[end]-finetruth(t[end])) < 1e-8
end

@testset "removed tolerance property cannot silently lose measured weights" begin
    probe = Parareal(RK4(h=0.01), RK4(h=0.1); criterion=WeightedDefect(weight=SecantWeight()))
    @test !hasproperty(probe, :tolerance)
    @test propertynames(probe) == fieldnames(typeof(probe))
    @test_throws (isdefined(Base, :FieldError) ? Base.FieldError : ErrorException) probe.tolerance
    @test_throws (isdefined(Base, :FieldError) ? Base.FieldError : ErrorException) (
        probe.tolerance.weights.w = 2.0
    )
    @test !occursin("tolerance", sprint(show, probe))
end
