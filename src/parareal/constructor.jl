"Immutable Parareal settings. Per-run propagation plans and weights live in the cache."
struct Parareal{
    F<:AbstractInitialValueSolver,
    G<:AbstractInitialValueSolver,
    P<:AbstractTimeParallelParameters,
    C<:AbstractCriterion,
} <: AbstractTimeParallelSolver
    finesolver::F
    coarsesolver::G
    parameters::P
    criterion::C
end
function Parareal(
    f::AbstractInitialValueSolver,
    g::AbstractInitialValueSolver;
    N::Integer=10,
    K::Integer=N,
    parameters=nothing,
    criterion::AbstractCriterion=RelativeUpdate(),
)
    p = parameters === nothing ? PararealParameters(; N, K) : parameters
    return Parareal(f, g, p, criterion)
end
"Build changed Parareal settings without changing the original."
function NSDEBase.remake(
    p::Parareal;
    finesolver=p.finesolver,
    coarsesolver=p.coarsesolver,
    criterion=p.criterion,
    N=p.parameters.N,
    K=(p.parameters.K == p.parameters.N ? N : p.parameters.K),
)
    parameters =
        N == p.parameters.N && K == p.parameters.K ? p.parameters : PararealParameters(N, K)
    return Parareal(finesolver, coarsesolver; parameters, criterion)
end
