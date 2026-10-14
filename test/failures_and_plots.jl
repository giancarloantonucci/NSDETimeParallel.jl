@testset "Chunk exceptions become nonfinite returns" begin
    p=IVP((u, t)->u .^ 2, [1.0], (0.0, 0.4))
    for backend in (SerialBackend(), ThreadsBackend())
        for (fine, coarse) in (
            (BackwardEuler(h=0.1, Mₙ=1, εᵣ=1e-14), Euler(h=0.05)),
            (RK4(h=0.01), BackwardEuler(h=0.1, Mₙ=1, εᵣ=1e-14)),
        )
            solver=Parareal(fine, coarse; N=4, criterion=NoStop())
            result=@test_logs (:warn, r"Chunk solve failed") match_mode=:any solve(
                p, solver; backend
            )
            @test result.diagnostics.status == :nonfinite
            @test result.errors === result.diagnostics.residuals
            @test all(u->all(isnan, u), boundaryvalues(result.lastiterate))
        end
        # The fine solver exhausts nonfinite adaptive retries; Euler on the
        # chunk starts still succeeds, so failure occurs inside a fine batch.
        q=IVP(
            (u, t) -> (t==0 || t==0.1 || t==0.2 || t==0.3) ? -u : fill(NaN, length(u)),
            [1.0],
            (0.0, 0.4),
        )
        solver=Parareal(DormandPrince54(h=0.01, Mₙ=2), Euler(h=0.1); N=4, criterion=NoStop())
        result=@test_logs (:warn, r"Chunk solve failed") match_mode=:any solve(q, solver; backend)
        @test result.diagnostics.status == :nonfinite
    end
end

@testset "Convergence recipe keeps positive finite residuals and batch indices" begin
    sol=solve(problem, make_parareal(N=4))
    resize!(sol.errors, 5)
    sol.errors .= [1e-3, 0.0, Inf, NaN, 1e-9]
    series=RecipesBase.apply_recipe(Dict{Symbol,Any}(), NSDEBase._Convergence(sol))
    @test series[1].plotattributes[:yscale] == :log10
    @test series[1].args == ([1, 5], [1e-3, 1e-9])
    phase=RecipesBase.apply_recipe(
        Dict{Symbol,Any}(:label=>"run"), NSDEBase._PhasePlot(sol.lastiterate)
    )
    @test phase[1].plotattributes[:label] == "run"
    @test all(s->s.plotattributes[:label]=="", phase[2:end])
end

struct DeliberateSolverFailure <: SolverFailure end

struct FailOnCall{S} <: NSDEBase.AbstractInitialValueSolver
    solver::S
    limit::Int
end
function NSDEBase.initialize_cache(p::NSDEBase.AbstractInitialValueProblem, s::FailOnCall)
    return (inner=initialize_cache(p, s.solver), calls=Ref(0))
end
function NSDEBase.initialize_solution(
    p::NSDEBase.AbstractInitialValueProblem, s::FailOnCall; kwargs...
)
    return initialize_solution(p, s.solver; kwargs...)
end
function NSDEBase.solve!(c, sol, p::NSDEBase.AbstractInitialValueProblem, s::FailOnCall)
    c.calls[] += 1
    c.calls[] == s.limit && throw(DeliberateSolverFailure())
    return solve!(c.inner, sol, p, s.solver)
end
@testset "Failure during coarse correction, then cache reuse" begin
    p=Logistic(0.3, (0.0, 1.0))
    solver=Parareal(RK4(h=0.01), FailOnCall(Euler(h=0.1), 4); N=4, criterion=NoStop())
    for backend in (SerialBackend(), ThreadsBackend())
        cache=initialize_cache(p, solver)
        failed=@test_logs (:warn, r"Chunk solve failed") solver(cache, p; backend)
        @test failed.diagnostics.status == :nonfinite
        @test failed.diagnostics.fine_batches == 1
        @test failed.diagnostics.sweeps == 1
        # A new solver rebuilds its caches and clears the failure flags.
        good=remake(solver; coarsesolver=Euler(h=0.1))
        result=good(cache, p; backend)
        @test result.diagnostics.status == :exact
        @test boundaryvalues(result.lastiterate) == reference(p, good)
    end
end

@testset "Programming errors propagate outside MPI" begin
    # N=1 avoids the coarse pass, so threaded errors originate in a fine task.
    bad=IVP((u, t)->throw(BoundsError(u, 99)), [1.0], (0.0, 0.1))
    par=Parareal(RK4(h=0.01), Euler(h=0.02); N=1, criterion=NoStop())
    @test_throws BoundsError solve(bad, par; backend=SerialBackend())
    err = try
        solve(bad, par; backend=ThreadsBackend())
        nothing
    catch e
        e
    end
    @test err isa CompositeException
    @test occursin("BoundsError", sprint(showerror, err))
    # The coarse pass follows the same policy.
    @test_throws BoundsError solve(bad, remake(par; N=2); backend=SerialBackend())
end
