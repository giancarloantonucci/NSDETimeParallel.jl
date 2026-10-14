# Run with an MPI launcher compatible with the MPI.jl environment, for example:
# mpiexecjl -n 2 julia --project examples/mpi.jl
using MPI
using NSDERungeKutta, NSDETimeParallel
MPI.Init()
try
    comm = MPI.COMM_WORLD
    N = MPI.Comm_size(comm)
    problem = Logistic(0.3, (0.0, 2.0))
    par = Parareal(RK4(; h=0.005), RK4(; h=0.05); N=N, criterion=NoStop())
    result = solve(problem, par; backend=MPIBackend())
    @assert result.diagnostics.status == :exact
    @assert length(result.errors) == N
    @assert result.errors === result.diagnostics.residuals
    @assert boundaryvalues(result.lastiterate) == reference(problem, par)
    if MPI.Comm_rank(comm) == 0
        println((status=result.diagnostics.status, batches=result.diagnostics.fine_batches))
    end
finally
    MPI.Finalize()
end
