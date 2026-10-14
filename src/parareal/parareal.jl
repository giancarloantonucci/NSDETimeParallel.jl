# src/parareal/parareal.jl
#
# THE Parareal algorithm — written once. Backends (serial, threads,
# distributed, MPI) supply execution primitives only; see backends.jl. Any fix
# made here reaches every backend by construction.

"""
    coarseinit!(cache, parareal; mpi=false)

runs the serial coarse pass over the chunk grid, filling `G` and — where
`makeGs` allows — seeding the chunk starting values `U`. Chunks whose
`makeGs[n]` is `false` keep whatever `U[n]` already holds (injected guesses,
e.g. from a moving-window driver).
"""
function coarseinit!(cache::PararealCache, parareal::Parareal; mpi::Bool=false)
    @↓ makeGs, U, T, F, G, chunkproblems, coarsecache, coarsesolution = cache
    @↓ coarsesolver = cache
    N = parareal.parameters.N
    # F[1]/G[1] mirror the window's starting value; on cache reuse (window
    # drivers) the start has moved, so refresh them every run:
    copyto!(F[1], U[1])
    copyto!(G[1], U[1])
    for n in 1:(N - 1)
        cache.chunk_failed[n] |=
            !solve_chunk_kernel!(coarsecache, coarsesolution, chunkproblems[n], coarsesolver; mpi)
        copyto!(G[n + 1], chunk_endpoint(coarsesolution, T[n + 1]))
        if makeGs[n + 1]
            copyto!(U[n + 1], G[n + 1])
        end
    end
    return cache
end

"""
    correct!(cache, parareal, k; mpi=false)

the serial Parareal correction sweep after fine batch `k`: re-run the coarse
solver from the freshest starts and update `U[n+1] = F[n+1] + (Gnew − Gold)`,
writing in place against the `U_` snapshot maintained by the iteration loop.

Two details make finite termination exact in floating point, not just in real
arithmetic:

- **Order of operations.** `F + (Gnew − Gold)` returns `F` bitwise whenever the
  two coarse values cancel bitwise; `(Gnew + F) − Gold` does not — it rounds
  `F` into `Gnew` first and loses low bits (with `Gnew = Gold = 1e16` and
  `F = 1` it returns `0`). Do not "simplify" this expression.
- **The newly exact interface is assigned, not computed.** After batch `k` the
  start `U[k]` is already exact, so `F[k+1]` is the exact value of `U[k+1]`
  and `U[k]` has not moved since `G[k+1]` was last computed: the coarse re-run
  at `n = k` would reproduce `G[k+1]` bitwise. Skip it and copy `F[k+1]`
  across — one coarse solve saved per sweep, and no arithmetic on the exact
  prefix at all.
"""
function correct!(cache::PararealCache, parareal::Parareal, k::Integer; mpi::Bool=false)
    @↓ U, T, F, G, chunkproblems, coarsecache, coarsesolution = cache
    @↓ coarsesolver = cache
    N = parareal.parameters.N
    k ≤ N - 1 || return cache
    copyto!(U[k + 1], F[k + 1]) # exact: fine output from an exact start
    for n in (k + 1):(N - 1)
        cache.chunk_failed[n] |=
            !solve_chunk_kernel!(coarsecache, coarsesolution, chunkproblems[n], coarsesolver; mpi)
        v = chunk_endpoint(coarsesolution, T[n + 1])
        @. U[n + 1] = F[n + 1] + (v - G[n + 1])
        copyto!(G[n + 1], v)
    end
    return cache
end

"""
    parareal!(cache, solution, problem, parareal, backend; saveiterates, directory, nocollect)

runs Parareal with the given [`AbstractPararealBackend`](@ref). Under MPI this
is SPMD: every rank calls it; `is_root` guards the serial parts.
"""
function parareal!(
    cache::PararealCache,
    solution::PararealSolution,
    problem::AbstractInitialValueProblem,
    parareal::Parareal,
    backend::AbstractPararealBackend;
    saveiterates::Bool=false,
    directory::Union{Nothing,String}=nothing,
    nocollect::Bool=false,
    record_interfaces::Bool=false,
    initial_weight=nothing,
)
    @↓ N, K = parareal.parameters
    criterion = parareal.criterion
    prepare_chunks!(cache, parareal)
    cache.weightstate = WeightState(weight_rule(parareal.criterion))
    cache.weightstate.w =
        initial_weight === nothing ? weight_rule(parareal.criterion).w : initial_weight
    isfinite(cache.weightstate.w) && cache.weightstate.w > 0 ||
        throw(ArgumentError("initial weight must be finite and positive"))
    weights = cache.weightstate
    @↓ errors = solution

    if backend isa MPIBackend
        MPI.Initialized() ||
            throw(ErrorException("`MPIBackend` needs `MPI.Init()` before solving."))
        MPI.Comm_size(MPI.COMM_WORLD) >= N ||
            throw(ArgumentError("MPIBackend requires at least N MPI ranks"))
        if saveiterates
            # All ranks must agree on the run directory; disk is used ONLY for
            # this opt-in history, never for the core hand-off.
            directory = MPI.bcast(
                directory === nothing ? mktempdir() : directory, 0, MPI.COMM_WORLD
            )
        end
    end

    solution.diagnostics = PararealDiagnostics(problem.u0, weights.w, errors)
    diag = solution.diagnostics
    if is_root(backend)
        # A reused solution was resized down to its previous run's k_final;
        # restore full capacity so `errors[k]` (and the iterate history) can
        # be written again — reuse must behave exactly like a fresh solution.
        length(errors) < K && resize!(errors, K)
        if saveiterates && !(solution.iterates isa Nothing) && length(solution.iterates) < K
            append!(
                solution.iterates,
                [PararealIterate(problem, parareal) for _ in (length(solution.iterates) + 1):K],
            )
        end
        coarseinit!(cache, parareal; mpi=backend isa MPIBackend)
        diag.coarse_chunk_calls = N-1
        diag.initial_starts = deepcopy(cache.U)
    end

    k_final = 0
    for k in 1:min(K, N)
        starts_finite =
            is_root(backend) ? !any(cache.chunk_failed) && all(u->all(isfinite, u), cache.U) : true
        starts_finite = sync_converged(backend, starts_finite)
        if !starts_finite
            diag.status = :nonfinite
            invalidate_chunks!(solution.lastiterate, cache.T, k:N)
            break
        end
        sync_starts!(backend, cache, k)
        fine_map!(backend, cache, solution, problem, parareal, k; saveiterates, directory)
        k_final = k
        diag.fine_batches += 1
        diag.fine_chunk_calls += N-k+1
        # Root needs current terminal outputs for the secant updater and saved histories.
        if (record_interfaces || needs_fine_history(weights.rule)) && backend isa MPIBackend
            collect_chunks!(backend, solution, cache; saveiterates=false)
        end
        local_finite = if backend isa MPIBackend
            n = MPI.Comm_rank(MPI.COMM_WORLD)+1
            n > N || (!cache.chunk_failed[n] && all(isfinite, solution.lastiterate[n].u[end]))
        else
            !any(cache.chunk_failed) && all(n->all(isfinite, solution.lastiterate[n].u[end]), 1:N)
        end
        finite =
            backend isa MPIBackend ? MPI.Allreduce(local_finite, &, MPI.COMM_WORLD) : local_finite
        converged = false
        if is_root(backend)
            if !finite
                errors[k] = Inf
                push!(diag.weights, weights.w)
                diag.status = :nonfinite
                converged = true
            else
                update_weight!(weights, cache, solution, k)
                errors[k] = residual(criterion, cache, k)
                push!(diag.weights, weights.w)
                record_interfaces && push!(
                    diag.interfaces, [deepcopy(cache.U); [copy(solution.lastiterate[N].u[end])]]
                )
                saveiterates && snapshot!(backend, solution, k)
                if k == N
                    diag.status = :exact
                    converged = true
                elseif accept(criterion, errors[k])
                    diag.status = :converged
                    converged = true
                elseif k == K
                    diag.status = :budget_exhausted
                    converged = true
                end
            end
        end
        converged = sync_converged(backend, converged)
        converged && break
        if is_root(backend)
            # Snapshot independently of ψ: weighting and update criteria share
            # the same old starts, and ψ₂ does not itself maintain U_.
            for n in 1:N
                copyto!(cache.U_[n], cache.U[n])
            end
            correct!(cache, parareal, k; mpi=backend isa MPIBackend)
            diag.sweeps += 1
            diag.coarse_chunk_calls += N-k-1
        end
    end
    if is_root(backend)
        resize!(errors, k_final)
        if saveiterates && !(solution.iterates isa Nothing)
            resize!(solution.iterates, k_final)
        end
    end

    if backend isa MPIBackend
        solution.diagnostics = MPI.bcast(is_root(backend) ? diag : nothing, 0, MPI.COMM_WORLD)
        solution.errors = solution.diagnostics.residuals
    end

    if !nocollect
        collect_chunks!(backend, solution, cache; saveiterates, directory)
    end
    return solution
end
