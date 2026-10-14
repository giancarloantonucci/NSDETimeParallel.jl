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
    exact = solve(problem, make_parareal(N=4, K=4, criterion=NoStop()))
    @test maxseam(exact) == 0.0
    early = solve(problem, make_parareal(N=4, K=1, criterion=NoStop()))     # one sweep: visible seams
    @test all(≥(0), seams(early))
    @test maxseam(early) > maxseam(full)
    single = solve(problem, make_parareal(N=1, K=1))
    @test isempty(seams(single))
    @test maxseam(single) == 0.0
end

@testset "flatten" begin
    N = 4
    solution = solve(problem, make_parareal(N=N))
    u, t = flatten(solution)
    @test length(u) == length(t)
    @test length(t) == sum(length(solution[n].t) for n in 1:N) - (N - 1) # boundaries deduplicated
    @test issorted(t)
    @test u[1] == problem.u0
    @test norm(u[end] - finetruth(t[end])) < 1e-8
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
    vecsolution = solve(
        vecproblem,
        Parareal(
            finesolver,
            coarsesolver;
            parameters=PararealParameters(N=3, K=3),
            criterion=MeanRelativeUpdate(ϵ=1e-10),
        ),
    )
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

@testset "RecipesBase recipes (headless)" begin
    psol = solve(problem, make_parareal(N=4))
    for obj in (
        psol,
        psol.lastiterate,
        NSDEBase._PhasePlot(psol),
        NSDEBase._PhasePlot(psol.lastiterate),
        NSDEBase._Convergence(psol),
    )
        @test RecipesBase.apply_recipe(Dict{Symbol,Any}(), obj) isa Vector{RecipesBase.RecipeData}
    end
end

@testset "Wnorm uses measured weights only" begin
    p = Lorenz([2.0, 3.0, -14.0], (0.0, 1.0))
    solver = Parareal(
        RK4(h=0.002), RK4(h=0.02); N=8, K=2, criterion=WeightedDefect(ϵ=0.0, weight=SecantWeight())
    )
    sol = solve(p, solver)
    ref = reference(p, solver)
    w = sol.diagnostics.weights[end]
    @test w > 1
    @test Wnorm(sol, ref) == Wnorm(sol.lastiterate, ref, w)
    @test Wnorm(sol, ref) == Wnorm(sol.lastiterate, ref, FixedWeight(w))
    @test Wnorm(sol, ref) < Wnorm(sol.lastiterate, ref, 1.0)
    @test_throws ArgumentError Wnorm(sol.lastiterate, ref, SecantWeight())
    @test_throws ArgumentError Wnorm(sol.lastiterate, ref, NeighbourWeight())
    @test make_parareal().criterion isa RelativeUpdate
end
