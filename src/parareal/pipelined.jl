"""
    PipelinedMPIBackend()

MPI frontier scheduling with one rank per chunk. Supported criteria are
`WeightedDefect(p=Inf, weight=FixedWeight(...))`, `FrontierRelativeUpdate`,
and `NoStop`. Other criteria are rejected; the bulk all-interface update
rule and changing secant weights require synchronized iterates.

Maximum-defect acceptance waits for the upstream start to freeze before
accepting its outgoing defect. Thus the assembled accepted starts satisfy
the selected defect tolerance. Converting that tolerance to a trajectory
error additionally requires bounds on the fine-map amplification.

Call `MPI.Init()` first and launch at least `N` ranks. Global histories are
unsupported. Diagnostics report actual total fine/coarse calls and per-chunk
`local_sweeps`; `sweeps` and `fine_batches` are `nothing` because this schedule
has no synchronized batches. `residuals` contains the maximum diagnostic
among ranks active at each local iteration index, not a bulk residual.

The final gather returns each fine trajectory from its accepted start.
With `NoStop()` and `K >= N`, completion reproduces `reference(problem, solver)`.
"""
struct PipelinedMPIBackend <: AbstractPararealBackend end

is_root(::PipelinedMPIBackend) = MPI.Comm_rank(MPI.COMM_WORLD) == 0

function parareal!(
    cache::PararealCache,
    solution::PararealSolution,
    problem::AbstractInitialValueProblem,
    parareal::Parareal,
    ::PipelinedMPIBackend;
    saveiterates::Bool=false,
    directory::Union{Nothing,String}=nothing,
    nocollect::Bool=false,
)
    saveiterates && throw(
        ArgumentError(
            "`PipelinedMPIBackend` does not support `saveiterates`: recording every sweep forces a per-sweep synchronisation — exactly the serial cost pipelining removes. Use `MPIBackend` for history runs.",
        ),
    )
    MPI.Initialized() ||
        throw(ErrorException("`PipelinedMPIBackend` needs `MPI.Init()` before solving."))

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    @↓ N, K = parareal.parameters
    ϵ = threshold(parareal.criterion)
    MPI.Comm_size(comm) ≥ N || throw(
        ErrorException(
            "`PipelinedMPIBackend` needs one rank per chunk: N = $N chunks but only $(MPI.Comm_size(comm)) ranks.",
        ),
    )
    kind = pipeline_kind(parareal.criterion)
    useψ∞ = kind === :defect
    solution.diagnostics = PararealDiagnostics(
        problem.u0, weight_rule(parareal.criterion).w, solution.errors
    )
    fine_calls = 0
    coarse_calls = 0

    @↓ U, T, chunkproblems, coarsecache, coarsesolution = cache
    @↓ coarsesolver = cache
    n = rank + 1            # the chunk this rank owns
    kdone = 0
    errhist = Float64[]

    if n ≤ N
        # ---------------- pipelined coarse initialisation ----------------
        # The chain replaces the root's serial coarseinit!: rank n receives
        # its start from the left, propagates one chunk, hands off — and is
        # then free to begin sweep 1's fine solve while the chain continues
        # downstream. Tags stay tiny (≤ 2K + 1): the MPI standard only
        # guarantees tags up to 32767.
        if n == 1
            copyto!(U[1], problem.u0)
        else
            MPI.Recv!(U[n], comm; source=rank - 1, tag=0)
        end
        cache.chunk_failed[n] |=
            !solve_chunk_kernel!(
                coarsecache, coarsesolution, chunkproblems[n], coarsesolver; mpi=true
            )
        coarse_calls += 1
        Gold = copy(chunk_endpoint(coarsesolution, T[n + 1]))       # G(Uₙ⁰) over the own chunk
        sendval = copy(Gold)
        flagsend = zeros(UInt8, 1)
        reqs = MPI.Request[]
        if n < N
            push!(reqs, MPI.Isend(sendval, comm; dest=rank + 1, tag=0))
        end

        # ------------------------- frontier loop -------------------------
        left_converged = (n == 1)   # chunk 1's start is exact from the outset
        converged = false
        uold = copy(Gold)  # the boundary value sent last sweep; at sweep 1 the coarse
        # prediction — which IS the downstream start U[n+1]⁰, so the
        # ψ∞ test at sweep 1 measures the true iterate-0 defect.
        unew = similar(Gold)
        flagrecv = zeros(UInt8, 1)
        dirty = false      # a start arrived after the last fine solve
        # ψ∞ mode: weight at THIS rank's outgoing boundary (start n + 1), per-unit-time
        # base as in ψ₂; rank N owns no Parareal unknown (the terminal boundary is
        # omitted by convention), so its own test is vacuous and it accepts on the flag.
        w∞ = useψ∞ ? max(1.0, weight_rule(parareal.criterion).w) : 1.0
        Wexp = (useψ∞ && n < N) ? w∞^(T[1] - T[n + 1]) : 0.0
        startprev = similar(U[n]) # detects a start CHANGING in the flag-raising message
        kmax = min(K, n)   # diagonal: chunk n's start is exact after n − 1 corrections
        k = 0
        while k < kmax && !converged
            k += 1
            # (1) Fine solve from the CURRENT start. This is the work that
            # hides the serial coarse chain: it runs BEFORE blocking on the
            # left neighbour — the entire pipelining gain lives here.
            solvechunk_local!(cache, solution, parareal, n)
            fine_calls += 1
            Fval = chunk_endpoint(solution.lastiterate[n], T[n + 1])
            # ψ∞ mode: the CURRENT defect of iterate k − 1 at this rank's
            # outgoing boundary — `uold` is the value sent last sweep, `Fval`
            # the fine solve from the very start that produced it: a
            # one-stage-lag evaluation on consistent starts.
            tψ = (useψ∞ && n < N) ? Wexp * norm(uold - Fval) : 0.0
            dirty = false
            startchanged = false
            # (2) Only now block on the corrected start and frontier flag:
            if !left_converged
                copyto!(startprev, U[n])
                MPI.Recv!(flagrecv, comm; source=rank - 1, tag=2k + 1)
                MPI.Recv!(U[n], comm; source=rank - 1, tag=2k)
                left_converged = !iszero(flagrecv[1])
                startchanged = U[n] != startprev
                dirty = true
            end
            # ψ∞ acceptance — decided BEFORE the correction. `tψ` tested
            # iterate k − 1, so acceptance is sound only if the start that
            # produced `uold`/`Fval` is the start now frozen: when the
            # flag-raising message also CHANGED the start (upstream diagonal
            # exit), defer one sweep and re-test against the frozen start.
            accepted = useψ∞ && left_converged && !startchanged && tψ ≤ ϵ
            if accepted
                push!(errhist, tψ)
                converged = true
                # The accepted iterate is k − 1: `lastiterate[n]` (the fine
                # solve from its start) is final, the accepted boundary is
                # `uold`, and this sweep's correction is discarded. `dirty`
                # is cleared even if a newer start arrived — re-solving from
                # it would break the certificate.
                dirty = false
            else
                # (3) Local coarse step from the (possibly updated) start,
                # then the Parareal correction — the RIGHT neighbour's next
                # start:
                cache.chunk_failed[n] |=
                    !solve_chunk_kernel!(
                        coarsecache, coarsesolution, chunkproblems[n], coarsesolver; mpi=true
                    )
                coarse_calls += 1
                Gnew = chunk_endpoint(coarsesolution, T[n + 1])
                @. unew = Fval + (Gnew - Gold) # parenthesised on purpose: see `correct!`
                copyto!(Gold, Gnew)
                # (4) Frontier convergence on the OUTGOING boundary:
                if useψ∞
                    push!(errhist, tψ)
                    converged = k == n # diagonal only; ψ∞ acceptance is above
                else
                    push!(errhist, norm(unew - uold) / max(norm(unew), eps()))
                    converged = (kind === :frontier && left_converged && errhist[end] ≤ ϵ) || k == n
                end
                copyto!(uold, unew)
            end
            # (5) Hand off. A `true` flag is, by construction, the LAST
            # message this rank ever sends; the receiver stops receiving on
            # consuming it, so send and receive counts match on every exit
            # path (acceptance, diagonal, K cap) and the protocol cannot
            # deadlock. On ψ∞ acceptance the value sent is the ACCEPTED
            # boundary `uold` — bitwise the value sent last sweep, so the
            # receiver's running fine solve stays valid; on the diagonal it
            # is the exact `unew`, and `startchanged` downstream defers
            # acceptance until the test has run against it.
            if n < N
                foreach(MPI.Wait, reqs)
                empty!(reqs)
                copyto!(sendval, accepted ? uold : unew)
                flagsend[1] = UInt8(converged)
                push!(reqs, MPI.Isend(flagsend, comm; dest=rank + 1, tag=2k + 1))
                push!(reqs, MPI.Isend(sendval, comm; dest=rank + 1, tag=2k))
            end
        end
        # The accepted chunk must be the fine solve FROM the accepted start:
        # if a start arrived after the last sweep's fine solve, run one more.
        # At the diagonal this makes the chunk the exact chunked fine solve.
        if dirty
            solvechunk_local!(cache, solution, parareal, n)
            fine_calls += 1
        end
        foreach(MPI.Wait, reqs)
        kdone = k
    end

    if n <= N && cache.chunk_failed[n]
        invalidate_chunks!(solution.lastiterate, T, n:n)
    end

    # -------------------- collection & diagnostics --------------------
    kfin = MPI.Allreduce(kdone, MPI.MAX, comm)
    histories = MPI.gather(errhist, comm; root=0)
    errors = MPI.bcast(
        if rank == 0
            [maximum((h[k] for h in histories if length(h) ≥ k); init=0.0) for k in 1:kfin]
        else
            nothing
        end,
        0,
        comm,
    )
    map!(x -> isfinite(x) ? x : Inf, errors, errors)
    resize!(solution.errors, kfin)
    copyto!(solution.errors, errors)
    diag = solution.diagnostics
    diag.fine_batches = nothing
    diag.sweeps = nothing
    diag.fine_chunk_calls = MPI.Allreduce(fine_calls, +, comm)
    diag.coarse_chunk_calls = MPI.Allreduce(coarse_calls, +, comm)
    diag.local_sweeps = MPI.Allgather(kdone, comm)[1:N]
    finite = MPI.Allreduce(
        n > N || (!cache.chunk_failed[n] && all(isfinite, solution.lastiterate[n].u[end])), &, comm
    )
    completed = MPI.Allreduce(n > N || converged, &, comm)
    diagonal = all(diag.local_sweeps[n] >= n for n in 1:N)
    diag.status = if !finite
        :nonfinite
    elseif completed
        (diagonal ? :exact : :converged)
    else
        :budget_exhausted
    end
    if !nocollect
        collect_chunks!(MPIBackend(), solution, cache)
    end
    return solution
end
