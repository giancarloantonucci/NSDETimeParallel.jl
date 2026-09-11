# MPI backend check — run OUTSIDE CI, on a machine with MPI:
#
#     julia --project -e 'using MPI; MPI.install_mpiexecjl()'   # once
#     mpiexecjl -n 4 julia --project=@. test/mpi/runtests_mpi.jl
#
# Every rank runs this SPMD; assertions execute on the root.
using MPI
MPI.Init()
using NSDETimeParallel, NSDERungeKutta, LinearAlgebra, Test

problem = Logistic(0.3, (0.0, 1.0))
finesolver = RK4(h=1e-3)
N = MPI.Comm_size(MPI.COMM_WORLD)
parareal = Parareal(finesolver, RungeKutta4(h=5e-2);
                    parameters=PararealParameters(N=N, K=N), tolerance=Tolerance(ϵ=1e-12))
solution = solve(problem, parareal; backend=MPIBackend())

if MPI.Comm_rank(MPI.COMM_WORLD) == 0
    reference = solve(problem, parareal)  # serial, same skeleton
    @testset "MPI ≡ SERIAL" begin
        @test solution.errors == reference.errors
        for n = 1:N
            @test solution[n].u == reference[n].u
        end
    end
end

# ------------------------- PipelinedMPIBackend -------------------------
# The pipelined backend is a SEMANTIC VARIANT (frontier convergence, ψ
# ignored), so "≡ SERIAL" is NOT the right assertion in general. What both
# algorithms share is the fixed point: at ϵ = 0, K = N, finite termination
# makes each produce the chunked fine solve, bitwise.
exact = Parareal(finesolver, RungeKutta4(h=5e-2);
                 parameters=PararealParameters(N=N, K=N), tolerance=Tolerance(ϵ=0.0))
psolution = solve(problem, exact; backend=PipelinedMPIBackend())
if MPI.Comm_rank(MPI.COMM_WORLD) == 0
    preference = solve(problem, exact)    # serial skeleton, same fixed point
    @testset "PIPELINED fixed point ≡ SERIAL (ϵ = 0, K = N)" begin
        @test numiterates(psolution) ≤ N  # frontier can only stop early
        for n = 1:N
            @test psolution[n].u == preference[n].u
            @test psolution[n].t == preference[n].t
        end
        @test maxseam(psolution) == maxseam(preference)
    end
end

# ---------------------- PipelinedMPIBackend + ψ∞ ----------------------
# The certified thesis criterion (§3.4): per-chunk weighted defects with
# prefix acceptance riding the pipeline. On acceptance the assembled vector
# must satisfy ψ∞ ≤ ϵ — asserted on the root by RECOMPUTING every weighted
# defect from the collected chunks with fresh serial fine solves. The weight
# must be preset; the updater is rejected outright.
ϵψ = 1e-8
ŵ = 1.0 # Logistic contracts; the unweighted base keeps the test strict everywhere
ψrun = Parareal(finesolver, RungeKutta4(h=5e-2);
                parameters=PararealParameters(N=N, K=N),
                tolerance=Tolerance(ϵ=ϵψ, ψ=ψ∞, weights=Weights(w=ŵ)))
ψsolution = solve(problem, ψrun; mode="PIPELINED")
if MPI.Comm_rank(MPI.COMM_WORLD) == 0
    @testset "PIPELINED honours ψ∞ (certified per-chunk defects)" begin
        Ts = boundarytimes(ψsolution.lastiterate)
        starts = [ψsolution[n].u[begin] for n = 1:N]
        for n = 1:N-1 # defect at start n+1: ‖U[n+1] − F(U[n])‖, weighted
            chunk = copy(problem, starts[n], Ts[n], Ts[n+1])
            Fend = solve(chunk, finesolver)(Ts[n+1])
            @test max(1.0, ŵ)^(Ts[1] - Ts[n+1]) * norm(starts[n+1] - Fend) ≤ ϵψ
        end
        @test numiterates(ψsolution) ≤ N
        @test all(e -> isfinite(e) && e ≥ 0, ψsolution.errors)
    end
    @testset "ψ∞ + updatew is rejected on the pipeline" begin
        bad = Parareal(finesolver, RungeKutta4(h=5e-2);
                       parameters=PararealParameters(N=N, K=N),
                       tolerance=Tolerance(ϵ=ϵψ, ψ=ψ∞, weights=Weights(w=ŵ, updatew=true)))
        @test_throws ArgumentError solve(problem, bad; mode="PIPELINED")
    end
end

# At a realistic tolerance the frontier must still land on the fine solve,
# and stop no later than the diagonal.
tolrun = Parareal(finesolver, RungeKutta4(h=5e-2);
                  parameters=PararealParameters(N=N, K=N), tolerance=Tolerance(ϵ=1e-10))
tsolution = solve(problem, tolrun; mode="PIPELINED") # the mode string dispatches too
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
MPI.Finalize()
