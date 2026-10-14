# MPI backend check — run collectively with MPI (also covered by CI):
#
#     julia --project -e 'using MPI; MPI.install_mpiexecjl()'   # once
#     mpiexecjl -n 4 julia --project=@. test/mpi/runtests_mpi.jl
#
# Every rank runs this SPMD; assertions execute on the root.
using MPI
MPI.Init()
using NSDETimeParallel, NSDERungeKutta, LinearAlgebra, Test

problem = Logistic(0.3, (0.0, 1.0))
finesolver = RK4(; h=1e-3)
N = MPI.Comm_size(MPI.COMM_WORLD)
parareal = Parareal(
    finesolver,
    RungeKutta4(; h=5e-2);
    parameters=PararealParameters(; N=N, K=N),
    criterion=RelativeUpdate(; atol=1e-12, rtol=1e-12),
)
solution = solve(problem, parareal; backend=MPIBackend())

if MPI.Comm_rank(MPI.COMM_WORLD) == 0
    serial_result = solve(problem, parareal)  # serial, same skeleton
    @testset "MPI ≡ SERIAL" begin
        @test solution.errors == serial_result.errors
        for n in 1:N
            @test solution[n].u == serial_result[n].u
        end
    end
end

# ------------------------- PipelinedMPIBackend -------------------------
# The pipelined backend is a SEMANTIC VARIANT (frontier convergence, ψ
# ignored), so "≡ SERIAL" is NOT the right assertion in general. What both
# algorithms share is the fixed point: at ϵ = 0, K = N, finite termination
# makes each produce the chunked fine solve, bitwise.
exact = Parareal(finesolver, RungeKutta4(; h=5e-2); N=N, K=N, criterion=NoStop())
psolution = solve(problem, exact; backend=PipelinedMPIBackend())
if MPI.Comm_rank(MPI.COMM_WORLD) == 0
    preference = solve(problem, exact)    # serial skeleton, same fixed point
    @testset "PIPELINED fixed point ≡ SERIAL (ϵ = 0, K = N)" begin
        @test numiterates(psolution) ≤ N  # frontier can only stop early
        for n in 1:N
            @test psolution[n].u == preference[n].u
            @test psolution[n].t == preference[n].t
        end
        @test maxseam(psolution) == maxseam(preference)
    end
end

# ---------------------- PipelinedMPIBackend + ψ∞ ----------------------
# The maximum weighted-defect criterion: per-chunk weighted defects with
# prefix acceptance riding the pipeline. On acceptance the assembled vector
# must satisfy ψ∞ ≤ ϵ — asserted on the root by RECOMPUTING every weighted
# defect from the collected chunks with fresh serial fine solves. The weight
# must be preset; the updater is rejected outright.
ϵψ = 1e-8
ŵ = 1.0 # Logistic contracts; the unweighted base keeps the test strict everywhere
ψrun = Parareal(
    finesolver,
    RungeKutta4(; h=5e-2);
    parameters=PararealParameters(; N=N, K=N),
    criterion=WeightedDefect(; ϵ=ϵψ, p=Inf, weight=FixedWeight(ŵ)),
)
ψsolution = solve(problem, ψrun; backend=PipelinedMPIBackend())
if MPI.Comm_rank(MPI.COMM_WORLD) == 0
    @testset "PIPELINED honours ψ∞ (certified per-chunk defects)" begin
        Ts = boundarytimes(ψsolution.lastiterate)
        starts = [ψsolution[n].u[begin] for n in 1:N]
        for n in 1:(N - 1) # defect at start n+1: ‖U[n+1] − F(U[n])‖, weighted
            chunk = remake(problem; u0=starts[n], tspan=(Ts[n], Ts[n + 1]))
            Fend = solve(chunk, finesolver)(Ts[n + 1])
            @test max(1.0, ŵ)^(Ts[1] - Ts[n + 1]) * norm(starts[n + 1] - Fend) ≤ ϵψ
        end
        @test numiterates(ψsolution) ≤ N
        @test all(e -> isfinite(e) && e ≥ 0, ψsolution.errors)
    end
    @testset "ψ∞ + adaptive weights is rejected on the pipeline" begin
        bad = Parareal(
            finesolver,
            RungeKutta4(h=5e-2);
            parameters=PararealParameters(N=N, K=N),
            criterion=WeightedDefect(ϵ=ϵψ, p=Inf, weight=NeighbourWeight(w=ŵ)),
        )
        @test_throws ArgumentError solve(problem, bad; backend=PipelinedMPIBackend())
    end
end

# At a realistic tolerance the frontier must still land on the fine solve,
# and stop no later than the diagonal.
tolrun = Parareal(
    finesolver, RungeKutta4(; h=5e-2); N=N, K=N, criterion=FrontierRelativeUpdate(; ϵ=1e-10)
)
tsolution = solve(problem, tolrun; backend=PipelinedMPIBackend())
if MPI.Comm_rank(MPI.COMM_WORLD) == 0
    finetruth = solve(problem, finesolver)
    @testset "PIPELINED converges to the fine solve" begin
        u, t = flatten(tsolution)
        @test issorted(t)
        @test norm(u[end] - finetruth(t[end])) < 1e-8
        @test numiterates(tsolution) ≤ N
        @test all(e -> isfinite(e) && e ≥ 0, tsolution.errors) # frontier diagnostic, not ψ
    end
end
# New bulk criteria, histories, and accounting must agree across backends.
for criterion in (:relative, :same_chunk)
    function make_tolerance()
        return if criterion == :relative
            RelativeUpdate(; atol=1e-12, rtol=1e-9)
        else
            WeightedDefect(; ϵ=1e-9, weight=SecantWeight())
        end
    end
    function make_solver()
        return Parareal(
            chunk_solver(RK4(; h=0.001), 1/N),
            chunk_solver(Euler(; h=0.01), 1/N);
            parameters=PararealParameters(; N=N),
            criterion=make_tolerance(),
        )
    end
    result=solve(problem, make_solver(); backend=MPIBackend(), history=:interfaces)
    if MPI.Comm_rank(MPI.COMM_WORLD)==0
        truth=solve(problem, make_solver(); history=:interfaces)
        @testset "Bulk $criterion across MPI and serial" begin
            @test result.errors==truth.errors
            @test result.diagnostics.status==truth.diagnostics.status
            @test result.diagnostics.sweeps==truth.diagnostics.sweeps
            @test result.diagnostics.fine_batches==truth.diagnostics.fine_batches
            @test result.diagnostics.coarse_chunk_calls==truth.diagnostics.coarse_chunk_calls
            @test result.diagnostics.interfaces==truth.diagnostics.interfaces
            @test boundaryvalues(result.lastiterate)==boundaryvalues(truth.lastiterate)
        end
    end
end

@testset "Pipeline diagnostics and capabilities" begin
    @test psolution.diagnostics.status == :exact
    @test psolution.diagnostics.fine_batches === nothing
    @test psolution.diagnostics.sweeps === nothing
    @test psolution.diagnostics.local_sweeps == collect(1:N)
    @test psolution.diagnostics.fine_chunk_calls >= N
    @test psolution.diagnostics.coarse_chunk_calls >= N
    @test psolution.diagnostics.residuals === psolution.errors
    unsupported = Parareal(finesolver, RK4(h=0.05); N, criterion=RelativeUpdate())
    @test_throws ArgumentError solve(problem, unsupported; backend=PipelinedMPIBackend())
end
@testset "Pipeline spare ranks and budget status" begin
    for chunks in (1, 2)
        p = Logistic(0.3, (2.0, 2.69))
        solver = Parareal(RK4(h=0.03), Euler(h=0.1); N=chunks, criterion=NoStop())
        result = solve(p, solver; backend=PipelinedMPIBackend())
        @test result.diagnostics.status == :exact
        @test result.diagnostics.local_sweeps == collect(1:chunks)
        if MPI.Comm_rank(MPI.COMM_WORLD) == 0
            @test boundaryvalues(result.lastiterate) == reference(p, solver)
        end
    end
    capped = Parareal(RK4(h=0.03), Euler(h=0.1); N=N, K=1, criterion=NoStop())
    result = solve(problem, capped; backend=PipelinedMPIBackend())
    @test result.diagnostics.status == :budget_exhausted
    @test result.diagnostics.local_sweeps == ones(Int, N)
end
@testset "MPI nonfinite coarse pass" begin
    p = Dahlquist([1.0], (0.0, 16.0); λ=-1000.0)
    solver = Parareal(RK4(h=0.001), Euler(h=0.02); N)
    result = solve(p, solver; backend=MPIBackend())
    @test result.diagnostics.status == :nonfinite
    @test result.diagnostics.fine_batches == 0
    if MPI.Comm_rank(MPI.COMM_WORLD) == 0
        @test isempty(result.errors)
        @test all(u -> all(isnan, u), boundaryvalues(result.lastiterate))
        @test all(isnan, result(7.5))
    end
end
@testset "Cost model uses bulk diagnostics and refuses pipeline counts" begin
    d = solution.diagnostics
    @test theoretical_speedup(d, 0.1) == N/(d.fine_batches + 0.1*d.coarse_chunk_calls)
    @test_throws ArgumentError theoretical_speedup(psolution, 0.1)
end
@testset "MPI validates rank count and reused problems" begin
    too_many = Parareal(RK4(h=0.01), Euler(h=0.1); N=N+1)
    @test_throws ArgumentError solve(problem, too_many; backend=MPIBackend())
    base = IVP(-1.0, [1.0], (0, 2))
    solver = Parareal(RK4(h=0.03), Euler(h=0.1); N, criterion=NoStop())
    for backend in (MPIBackend(), PipelinedMPIBackend())
        cache = NSDEBase.initialize_cache(base, solver)
        for q in (remake(base; u0=[0.7]), IVP(-2.0, [0.4], (0, 2)))
            result = solver(cache, q; backend)
            if MPI.Comm_rank(MPI.COMM_WORLD) == 0
                @test boundaryvalues(result.lastiterate) == reference(q, solver)
            end
        end
    end
    nanproblem = IVP((u, t)->fill(NaN, length(u)), [1.0], (0.0, 1.0))
    result = solve(nanproblem, Parareal(RK4(h=0.1), Euler(h=0.1); N=1); backend=MPIBackend())
    @test result.diagnostics.status == :nonfinite
    @test length(result.diagnostics.weights) == length(result.diagnostics.residuals) == 1
end
@testset "Bulk MPI returns trimmed shared residual storage on every rank" begin
    p = Logistic(0.3, (0.0, 1.0))
    solver = Parareal(RK4(h=0.01), Euler(h=0.1); N=N, criterion=WeightedDefect(ϵ=1e6))
    cache = initialize_cache(p, solver)
    result = initialize_solution(p, solver)
    solver(cache, result, p; backend=MPIBackend())
    @test result.diagnostics.fine_batches == 1
    @test length(result.errors) == length(result.diagnostics.residuals) == 1 < N
    @test result.errors === result.diagnostics.residuals
    @test result.errors == MPI.bcast(result.errors, 0, MPI.COMM_WORLD)
    @test all(isfinite, result.errors)
    solver = remake(solver; criterion=RelativeUpdate(atol=1e-12, rtol=1e-12))
    solver(cache, result, p; backend=MPIBackend())
    @test result.errors === result.diagnostics.residuals
    @test length(result.errors) == result.diagnostics.fine_batches
    @test result.errors == MPI.bcast(result.errors, 0, MPI.COMM_WORLD)
    bad = Dahlquist([1.0], (0.0, 16.0); λ=-1000.0)
    failed = solve(bad, Parareal(RK4(h=0.001), Euler(h=0.02); N); backend=MPIBackend())
    @test failed.diagnostics.status == :nonfinite
    @test failed.errors === failed.diagnostics.residuals
    @test isempty(failed.errors)
end
@testset "Rank-local solver exceptions do not strand MPI peers" begin
    p=IVP((u, t)->u .^ 2, [1.0], (0.0, 0.4))
    for backend in (MPIBackend(), PipelinedMPIBackend())
        for (fine, coarse) in (
            (BackwardEuler(h=0.1, Mₙ=1, εᵣ=1e-14), Euler(h=0.05)),
            (RK4(h=0.01), BackwardEuler(h=0.1, Mₙ=1, εᵣ=1e-14)),
        )
            result=solve(p, Parareal(fine, coarse; N=4, criterion=NoStop()); backend)
            @test result.diagnostics.status == :nonfinite
            @test result.errors === result.diagnostics.residuals
            @test all(x->isfinite(x)||isinf(x), result.errors)
            @test result.diagnostics.status ==
                MPI.bcast(result.diagnostics.status, 0, MPI.COMM_WORLD)
        end
        localbad=IVP(
            (u, t) -> (0.21<t<0.29 ? error("rank-local RHS failure") : -u), [1.0], (0.0, 0.4)
        )
        failed=solve(
            localbad, Parareal(RK4(h=0.01), Euler(h=0.1); N=4, criterion=NoStop()); backend
        )
        @test failed.diagnostics.status == :nonfinite
        @test failed.diagnostics.status == MPI.bcast(failed.diagnostics.status, 0, MPI.COMM_WORLD)
        q=IVP((u, t) -> (isinteger(10t) ? -u : fill(NaN, length(u))), [1.0], (0.0, 0.4))
        result=solve(
            q,
            Parareal(DormandPrince54(h=0.01, Mₙ=2), Euler(h=0.1); N=4, criterion=NoStop());
            backend,
        )
        @test result.diagnostics.status == :nonfinite
    end
end
MPI.Finalize()
