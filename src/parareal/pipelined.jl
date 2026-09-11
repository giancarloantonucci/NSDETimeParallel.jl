# src/parareal/pipelined.jl
#
# THE SEMANTIC VARIANT: pipelined (Aubanel-scheduled) Parareal over MPI with
# FRONTIER convergence. This is deliberately a SECOND algorithm, not a sixth
# adapter: pipelining fuses coarse propagation, fine solve and correction into
# one per-rank loop, which cannot be expressed with the five root-centric
# primitives of backends.jl. Everything that CAN be shared is shared — the
# cache, the chunk kernel `solvechunk_local!`, the final gather — and the loop
# below is the one thing that must differ. The semantic differences are stated
# loudly in the docstring; do not paper over them.

"""
    PipelinedMPIBackend <: AbstractPararealBackend

Pipelined (task-scheduled) Parareal over MPI ranks, one rank per chunk — the
scheduling of Aubanel (thesis §2.7, `eq:aubanel_parallel_efficiency`). There
is no root: the serial coarse chain is passed rank to rank, and each rank
runs its fine solve for the next sweep BEFORE blocking on the corrected start
from its left neighbour, so the fine work hides the serial coarse cost. Judge
measured speed-ups against
[`theoretical_speedup`](@ref)`(K, N, ζ; estimate = :aubanel)`.

This is a SEMANTIC VARIANT of Parareal, not a drop-in adapter. Differences
from the primitive-based backends:

- **Frontier convergence, two flavours.** Convergence sweeps left to right
  as a frontier; a chunk stops once every chunk to its left has accepted AND
  its own per-chunk test passes, or it hits the finite-termination diagonal
  (`k = n`). The per-chunk test depends on `tolerance.ψ`:

  * `ψ = ψ∞` — the CERTIFIED thesis criterion (§3.4 `subsec:local_proximity`):
    the test is the weighted CURRENT defect of the previous iterate at this
    rank's outgoing boundary, `w^(T[1] − T[n+1]) ‖uold − F(start)‖ ≤ ϵ`,
    where `uold` is the boundary value sent last sweep and `F(start)` this
    sweep's fine solve from the very start that produced it — the thesis's
    one-stage lag (`rem:pipelined_protocol`). The acceptance flag is the
    prefix conjunction `A_n = A_{n−1} ∧ tₙ`, riding the existing messages;
    when the run accepts, the assembled vector satisfies `ψ∞ ≤ ϵ`, hence
    `‖U − U*‖_{W,∞} ≤ s(θ_F) ϵ` (`thm:equivalence_error_norm_local`) — set
    `ϵ = r/s(θ_F)` to certify radius `r`. Requires a PRESET scalar weight:
    `weights.updatew = true` is rejected (no root sees all boundaries at a
    common sweep). Acceptance is suppressed for one sweep whenever the
    incoming start changed in the same message that raised the flag (an
    upstream diagonal exit): the test must be re-run against the frozen
    start, or the certificate would mix iterates — the thesis's freezing
    cascade, at most one extra sweep per chunk.
  * anything else — relative STAGNATION of the outgoing boundary
    (update ≤ `tolerance.ϵ`), the pre-thesis heuristic; the global `ψ` and
    its weights are ignored, with a warning unless `ψ = ψ₁`. Kept as the
    default for backward compatibility.

- **`solution.errors[k]`** records, among chunks still active at sweep `k`,
  the largest weighted defect (`ψ∞` mode: the running restriction of `ψ∞` to
  active chunks) or the largest relative update (stagnation mode) — a
  frontier diagnostic either way.
- **Chunks may stop at different sweeps**, so `numiterates(solution)` counts
  the sweeps of the LAST chunk to converge.
- **`saveiterates` is not supported**: recording every sweep would force the
  per-sweep synchronisation the pipeline exists to avoid.

Accepted chunks remain fine solves from their accepted starts (one final fine
solve is issued whenever a start arrived after the last sweep's), so
[`flatten`](@ref), [`seams`](@ref) and every downstream statistic read
exactly as for the other backends. At `ϵ = 0`, `K = N`, finite termination
makes the result the chunked fine solve, bitwise — the same fixed point as
`SerialBackend`, reached by a different route.

Requires `MPI.Init()` and at least `N` ranks; ranks beyond `N` idle through
the chain and join the final gather. Select with `mode = "PIPELINED"` or by
passing the backend object; pair with `Tolerance(ϵ = r/s, ψ = ψ∞,
weights = Weights(w = ŵ))` for the certified criterion (`ŵ` preset — the
updater is rejected here).
"""
struct PipelinedMPIBackend <: AbstractPararealBackend end

is_root(::PipelinedMPIBackend) = MPI.Comm_rank(MPI.COMM_WORLD) == 0

function parareal!(cache::PararealCache, solution::PararealSolution,
                   problem::AbstractInitialValueProblem, parareal::Parareal,
                   ::PipelinedMPIBackend;
                   saveiterates::Bool=false, directory::Union{Nothing,String}=nothing,
                   nocollect::Bool=false)
    saveiterates && throw(ArgumentError("`PipelinedMPIBackend` does not support `saveiterates`: recording every sweep forces a per-sweep synchronisation — exactly the serial cost pipelining removes. Use `MPIBackend` for history runs."))
    MPI.Initialized() || throw(ErrorException("`PipelinedMPIBackend` needs `MPI.Init()` before solving."))

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    @↓ N, K = parareal.parameters
    @↓ ϵ = parareal.tolerance
    MPI.Comm_size(comm) ≥ N ||
        throw(ErrorException("`PipelinedMPIBackend` needs one rank per chunk: N = $N chunks but only $(MPI.Comm_size(comm)) ranks."))
    useψ∞ = parareal.tolerance.ψ === ψ∞
    if useψ∞
        parareal.tolerance.weights.updatew &&
            throw(ArgumentError("`PipelinedMPIBackend` with `ψ∞` needs a PRESET weight: `update!` requires a root that sees all boundaries at a common sweep, which the pipeline abolishes. Measure `w` on a probe window (see the criterion docs) and pass `Weights(w = ŵ)`."))
        parareal.tolerance.weights.w isa Real ||
            throw(ArgumentError("`PipelinedMPIBackend` with `ψ∞` needs a SCALAR weight (the per-unit-time base `w`)."))
    elseif rank == 0 && parareal.tolerance.ψ !== ψ₁
        @warn "PipelinedMPIBackend with ψ ∉ (ψ₁, ψ∞) converges by a PER-BOUNDARY frontier (relative stagnation ≤ ϵ); this global ψ and its weights are ignored — a global criterion would re-serialise the pipeline every sweep. Use ψ∞ for the certified per-chunk criterion." maxlog = 1
    end

    @↓ U, T, chunkproblems, coarsecache, coarsesolution = cache
    @↓ coarsesolver = parareal
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
        NSDEBase.solve!(coarsecache, coarsesolution, chunkproblems[n], coarsesolver)
        Gold = copy(coarsesolution(T[n+1]))       # G(Uₙ⁰) over the own chunk
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
        w∞ = useψ∞ ? max(1.0, parareal.tolerance.weights.w) : 1.0
        Wexp = (useψ∞ && n < N) ? w∞^(T[1] - T[n+1]) : 0.0
        startprev = similar(U[n]) # detects a start CHANGING in the flag-raising message
        kmax = min(K, n)   # diagonal: chunk n's start is exact after n − 1 corrections
        k = 0
        while k < kmax && !converged
            k += 1
            # (1) Fine solve from the CURRENT start. This is the work that
            # hides the serial coarse chain: it runs BEFORE blocking on the
            # left neighbour — the entire pipelining gain lives here.
            solvechunk_local!(cache, solution, parareal, n)
            Fval = solution.lastiterate[n](T[n+1])
            # ψ∞ mode: the CURRENT defect of iterate k − 1 at this rank's
            # outgoing boundary — `uold` is the value sent last sweep, `Fval`
            # the fine solve from the very start that produced it: the
            # thesis's one-stage-lag evaluation (rem:pipelined_protocol).
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
                NSDEBase.solve!(coarsecache, coarsesolution, chunkproblems[n], coarsesolver)
                Gnew = coarsesolution(T[n+1])
                @. unew = Fval + (Gnew - Gold) # parenthesised on purpose: see `correct!`
                copyto!(Gold, Gnew)
                # (4) Frontier convergence on the OUTGOING boundary:
                if useψ∞
                    push!(errhist, tψ)
                    converged = k == n # diagonal only; ψ∞ acceptance is above
                else
                    push!(errhist, norm(unew - uold) / max(norm(unew), eps()))
                    converged = (left_converged && errhist[end] ≤ ϵ) || k == n
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
        dirty && solvechunk_local!(cache, solution, parareal, n)
        foreach(MPI.Wait, reqs)
        kdone = k
    end

    # -------------------- collection & diagnostics --------------------
    kfin = MPI.Allreduce(kdone, MPI.MAX, comm)
    histories = MPI.gather(errhist, comm; root=0)
    errors = MPI.bcast(
        rank == 0 ?
            [maximum((h[k] for h in histories if length(h) ≥ k); init=0.0) for k = 1:kfin] :
            nothing,
        0, comm)
    resize!(solution.errors, kfin)
    copyto!(solution.errors, errors)
    if !nocollect
        collect_chunks!(MPIBackend(), solution, cache)
    end
    return solution
end
