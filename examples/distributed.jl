using Distributed
project = dirname(Base.active_project())
workers_added = addprocs(2; exeflags="--project=$project")
try
    @everywhere using NSDERungeKutta, NSDETimeParallel
    problem = Logistic(0.3, (0.0, 2.0))
    par = Parareal(RK4(; h=0.005), RK4(; h=0.05); N=4, criterion=NoStop())
    result = solve(problem, par; backend=DistributedBackend())
    @assert boundaryvalues(result.lastiterate) == reference(problem, par)
    @assert result.diagnostics.status == :exact
finally
    rmprocs(workers_added)
end
