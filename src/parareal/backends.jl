"""
    AbstractPararealBackend

How the fine sweep of [`Parareal`](@ref) executes. The algorithm lives ONCE in
`parareal/parareal.jl`; a backend supplies only these primitives:

- `is_root(backend)` — whether this process runs the serial parts (coarse
  init, correction, convergence check). `true` everywhere except non-root MPI
  ranks.
- `sync_starts!(backend, cache, k)` — make the chunk starting values `U[n]`
  available wherever chunk `n` will be solved. No-op off MPI.
- `fine_map!(backend, cache, solution, problem, parareal, k; …)` — solve
  chunks `k:N` with the fine solver; deposit boundary values into
  `cache.F[n+1]` on the root and chunk solutions into `solution.lastiterate`
  (MPI ranks keep their chunk local until `collect_chunks!`).
- `snapshot!(backend, solution, k; …)` — record iterate `k` when
  `saveiterates` is on. In memory by default; on MPI the ranks have already
  written their chunk to the run directory inside `fine_map!`.
- `sync_converged(backend, flag)` — agree the convergence decision.
- `collect_chunks!(backend, solution, cache; …)` — bring final chunks (and,
  under MPI with `saveiterates`, the per-iterate history) to the root.
"""
abstract type AbstractPararealBackend end

"Runs everything in the calling task. The reference implementation."
struct SerialBackend <: AbstractPararealBackend end

"""
Fine chunks via `Threads.@threads` with dynamic scheduling in shared memory.
"""
struct ThreadsBackend <: AbstractPararealBackend end

"Fine chunks via `pmap` over the default worker pool (Distributed.jl)."
struct DistributedBackend <: AbstractPararealBackend end

"SPMD over MPI ranks: rank `r` owns chunk `r + 1`. Requires `MPI.Init()`."
struct MPIBackend <: AbstractPararealBackend end

# ------------------------------------------------------------------ defaults --

is_root(::AbstractPararealBackend) = true
sync_starts!(::AbstractPararealBackend, cache, k) = nothing
sync_converged(::AbstractPararealBackend, flag::Bool) = flag
collect_chunks!(::AbstractPararealBackend, solution, cache; kwargs...) = solution

function snapshot!(::AbstractPararealBackend, solution::PararealSolution, k::Integer; kwargs...)
    iterates = solution.iterates
    iterates isa Nothing && return solution
    for n in 1:numchunks(solution)
        iterates[k][n] = deepcopy(solution.lastiterate[n])
    end
    return solution
end

# The in-process chunk kernel shared by the serial and threaded backends: the
# innermost hot call. Everything it touches is preallocated; `solve!` reuses
# the chunk's solver cache and writes into the iterate's own chunk solution.
function solvechunk!(
    cache::PararealCache, solution::PararealSolution, parareal::Parareal, n::Integer
)
    @↓ T, F, chunkproblems, finecaches = cache
    cache.chunk_failed[n] |=
        !solve_chunk_kernel!(
            finecaches[n], solution.lastiterate[n], chunkproblems[n], cache.finesolver
        )
    if n < length(chunkproblems)
        copyto!(F[n + 1], chunk_endpoint(solution.lastiterate[n], T[n + 1]))
    end
    return nothing
end

# ------------------------------------------------------------------- serial --

function fine_map!(::SerialBackend, cache, solution, problem, parareal, k; kwargs...)
    for n in k:parareal.parameters.N
        solvechunk!(cache, solution, parareal, n)
    end
    return nothing
end

# ------------------------------------------------------------------ threads --

function fine_map!(::ThreadsBackend, cache, solution, problem, parareal, k; kwargs...)
    Threads.@threads for n in k:parareal.parameters.N # dynamic on ≥ 1.8 by default; see ThreadsBackend
        solvechunk!(cache, solution, parareal, n) # disjoint slots: no shared writes
    end
    return nothing
end

# -------------------------------------------------------------- distributed --

function fine_map!(::DistributedBackend, cache, solution, problem, parareal, k; kwargs...)
    @↓ U, T, F = cache
    @↓ finesolver = cache
    N = parareal.parameters.N
    # `pmap` over the DEFAULT pool: the schedule adapts to however many
    # workers exist, so the result cannot depend on the worker count. Chunk
    # solutions come back through memory (pmap return values), not disk.
    results = pmap(k:N) do n
        chunkproblem = remake(problem; u0=U[n], tspan=(T[n], T[n + 1]))
        chunksolution = NSDEBase.initialize_solution(chunkproblem, finesolver)
        work = NSDEBase.initialize_cache(chunkproblem, finesolver)
        ok = solve_chunk_kernel!(work, chunksolution, chunkproblem, finesolver)
        return (
            n=n, boundary=chunk_endpoint(chunksolution, T[n + 1]), chunk=chunksolution, failed=(!ok)
        )
    end
    for r in results
        solution.lastiterate[r.n] = r.chunk
        cache.chunk_failed[r.n] |= r.failed
        if r.n < N
            copyto!(F[r.n + 1], r.boundary)
        end
    end
    return nothing
end

# ---------------------------------------------------------------------- MPI --

is_root(::MPIBackend) = MPI.Comm_rank(MPI.COMM_WORLD) == 0

function sync_starts!(::MPIBackend, cache, k)
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    root = 0
    N = length(cache.U)
    if rank == root
        for dest in 1:(MPI.Comm_size(comm) - 1)
            w = dest + 1 # chunk owned by rank `dest`
            if k ≤ w ≤ N
                MPI.Send(cache.U[w], comm; dest, tag=w)
            end
        end
    else
        n = rank + 1
        if k ≤ n ≤ N
            MPI.Recv!(cache.U[n], comm; source=root, tag=n)
        end
    end
    return nothing
end

function fine_map!(
    ::MPIBackend,
    cache,
    solution,
    problem,
    parareal,
    k;
    saveiterates::Bool=false,
    directory::Union{Nothing,String}=nothing,
    kwargs...,
)
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    root = 0
    @↓ T, F = cache
    N = parareal.parameters.N
    n = rank + 1
    if k ≤ n ≤ N
        solvechunk_local!(cache, solution, parareal, n) # solve WITHOUT touching F (root-only state)
        if saveiterates && directory !== nothing
            filename = joinpath(directory, "iter_$(k)_chunk_$(n).jls")
            open(file -> serialize(file, (chunk_n=solution.lastiterate[n],)), filename, "w")
        end
        if n < N
            boundary = chunk_endpoint(solution.lastiterate[n], T[n + 1])
            if rank == root
                copyto!(F[n + 1], boundary)
            else
                MPI.Send(boundary, comm; dest=root, tag=N + n)
            end
        end
    end
    if rank == root
        for src in 1:(MPI.Comm_size(comm) - 1)
            w = src + 1
            if k ≤ w < N
                MPI.Recv!(F[w + 1], comm; source=src, tag=N + w)
            end
        end
    end
    return nothing
end

# Like `solvechunk!` but leaves `F` alone: under MPI only the root owns `F`.
function solvechunk_local!(
    cache::PararealCache, solution::PararealSolution, parareal::Parareal, n::Integer
)
    cache.chunk_failed[n] |=
        !solve_chunk_kernel!(
            cache.finecaches[n],
            solution.lastiterate[n],
            cache.chunkproblems[n],
            cache.finesolver;
            mpi=true,
        )
    return nothing
end

# Ranks already wrote their per-iterate chunks to disk inside `fine_map!`.
snapshot!(::MPIBackend, solution::PararealSolution, k::Integer; kwargs...) = solution

sync_converged(::MPIBackend, flag::Bool) = MPI.bcast(flag, 0, MPI.COMM_WORLD)

function collect_chunks!(
    ::MPIBackend,
    solution,
    cache;
    saveiterates::Bool=false,
    directory::Union{Nothing,String}=nothing,
    kwargs...,
)
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    root = 0
    N = numchunks(solution)
    n = rank + 1
    # Final chunk hand-off through memory (generic gather serialises objects),
    # not through the file system.
    payload = n ≤ N ? solution.lastiterate[n] : nothing
    gathered = MPI.gather(payload, comm; root)
    if rank == root
        for (r, chunk) in enumerate(gathered)
            m = r # gather order = rank order; rank r-1 owns chunk r
            if chunk !== nothing && m ≤ N
                solution.lastiterate[m] = chunk
            end
        end
        if saveiterates && directory !== nothing && !(solution.iterates isa Nothing)
            collect_iterates!(solution; directory)
        end
    end
    return solution
end

# Dynamic cache fields permit method/RHS changes. Specialize beyond this boundary.
@noinline function solve_chunk_kernel!(
    work::C, solution::S, problem::P, solver::F; mpi::Bool=false
) where {C,S,P,F}
    try
        NSDEBase.solve!(work, solution, problem, solver)
        return true
    catch err
        # Only numerical failures become nonfinite outside MPI. MPI must finish
        # its communication even after a rank-local programming error.
        if err isa InterruptException || err isa OutOfMemoryError || err isa StackOverflowError
            mpi &&
                MPI.Initialized() &&
                MPI.Comm_size(MPI.COMM_WORLD) > 1 &&
                MPI.Abort(MPI.COMM_WORLD, 1)
            rethrow()
        end
        (err isa SolverFailure || mpi) || rethrow()
        @warn "Chunk solve failed; returning :nonfinite" tspan=problem.tspan exception=(
            err, catch_backtrace()
        )
        invalidate_chunk!(solution, problem.tspan)
        return false
    end
end
