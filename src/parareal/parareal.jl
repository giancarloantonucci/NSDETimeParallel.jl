# src/parareal/parareal.jl
#
# THE Parareal algorithm — written once. Backends (serial, threads,
# distributed, MPI) supply execution primitives only; see backends.jl. Any fix
# made here reaches every backend by construction.

"""
    coarseinit!(cache, parareal)

runs the serial coarse pass over the chunk grid, filling `G` and — where
`makeGs` allows — seeding the chunk starting values `U`. Chunks whose
`makeGs[n]` is `false` keep whatever `U[n]` already holds (injected guesses,
e.g. from a moving-window driver).
"""
function coarseinit!(cache::PararealCache, parareal::Parareal)
    @↓ makeGs, U, T, F, G, chunkproblems, coarsecache, coarsesolution = cache
    @↓ coarsesolver = parareal
    N = parareal.parameters.N
    # F[1]/G[1] mirror the window's starting value; on cache reuse (window
    # drivers) the start has moved, so refresh them every run:
    copyto!(F[1], U[1])
    copyto!(G[1], U[1])
    for n = 1:N-1
        NSDEBase.solve!(coarsecache, coarsesolution, chunkproblems[n], coarsesolver)
        copyto!(G[n+1], coarsesolution(T[n+1]))
        if makeGs[n+1]
            copyto!(U[n+1], G[n+1])
        end
    end
    return cache
end

"""
    correct!(cache, parareal, k)

the serial Parareal correction sweep after fine batch `k`: re-run the coarse
solver from the freshest starts and update `U[n+1] = F[n+1] + (Gnew − Gold)`,
writing IN PLACE (ping-pong against the `U_` snapshot taken by the error
function — the old code rebound a fresh array here every step to protect that
history).

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
function correct!(cache::PararealCache, parareal::Parareal, k::Integer)
    @↓ U, T, F, G, chunkproblems, coarsecache, coarsesolution = cache
    @↓ coarsesolver = parareal
    N = parareal.parameters.N
    k ≤ N - 1 || return cache
    copyto!(U[k+1], F[k+1]) # exact: fine output from an exact start
    for n = k+1:N-1
        NSDEBase.solve!(coarsecache, coarsesolution, chunkproblems[n], coarsesolver)
        v = coarsesolution(T[n+1])
        @. U[n+1] = F[n+1] + (v - G[n+1])
        copyto!(G[n+1], v)
    end
    return cache
end

"""
    parareal!(cache, solution, problem, parareal, backend; saveiterates, directory, nocollect)

runs Parareal with the given [`AbstractPararealBackend`](@ref). Under MPI this
is SPMD: every rank calls it; `is_root` guards the serial parts.
"""
function parareal!(cache::PararealCache, solution::PararealSolution,
                   problem::AbstractInitialValueProblem, parareal::Parareal,
                   backend::AbstractPararealBackend;
                   saveiterates::Bool=false, directory::Union{Nothing,String}=nothing,
                   nocollect::Bool=false)
    @↓ N, K = parareal.parameters
    @↓ weights, ψ, ϵ = parareal.tolerance
    @↓ errors = solution

    if backend isa MPIBackend
        MPI.Initialized() || throw(ErrorException("`MPIBackend` needs `MPI.Init()` before solving."))
        if saveiterates
            # All ranks must agree on the run directory; disk is used ONLY for
            # this opt-in history, never for the core hand-off.
            directory = MPI.bcast(directory === nothing ? mktempdir() : directory, 0, MPI.COMM_WORLD)
        end
    end

    if is_root(backend)
        # A reused solution was resized down to its previous run's k_final;
        # restore full capacity so `errors[k]` (and the iterate history) can
        # be written again — reuse must behave exactly like a fresh solution.
        length(errors) < K && resize!(errors, K)
        if saveiterates && !(solution.iterates isa Nothing) && length(solution.iterates) < K
            append!(solution.iterates, [PararealIterate(problem, parareal) for _ = length(solution.iterates)+1:K])
        end
        coarseinit!(cache, parareal)
    end

    converged = false
    k_final = K
    for k = 1:K
        sync_starts!(backend, cache, k)
        fine_map!(backend, cache, solution, problem, parareal, k; saveiterates, directory)

        if is_root(backend)
            if saveiterates
                snapshot!(backend, solution, k)
            end
            update!(weights, cache.U, cache.F, cache.T)
            errors[k] = ψ(cache, k, weights)
            converged = errors[k] ≤ ϵ
        end
        converged = sync_converged(backend, converged)
        if converged
            k_final = k
            break
        end
        if is_root(backend)
            correct!(cache, parareal, k)
        end
    end

    if is_root(backend)
        resize!(errors, k_final)
        if saveiterates && !(solution.iterates isa Nothing)
            resize!(solution.iterates, k_final)
        end
    end

    if !nocollect
        collect_chunks!(backend, solution, cache; saveiterates, directory)
    end
    return solution
end
