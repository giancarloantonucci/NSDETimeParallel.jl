# Parareal scaling core: ONE measurement at this session's thread count.
# Driven by scaling_sweep.jl, which launches it at 1, 2, 4, 8, … threads and
# collates the curve. Can also be run directly:
#   julia --project=benchmark --threads=8 benchmark/scaling.jl [N] [K]
#
# What it measures
# ----------------
# Parareal wall time ≈ (coarse, SERIAL) + K × (fine chunk, PARALLEL over the
# pool) + K × (correction, SERIAL). To isolate the POOL we FIX the work: ϵ = 0
# and a fixed K, so every run does exactly K sweeps and the only variable is
# how many of the N fine chunks run at once. Only the fine sweep is parallel,
# so Amdahl caps the speed-up — seeing that ceiling is the point.
using NSDETimeParallel, NSDERungeKutta, Printf
using LinearAlgebra

const N = length(ARGS) ≥ 1 ? parse(Int, ARGS[1]) : 16
const K = length(ARGS) ≥ 2 ? parse(Int, ARGS[2]) : 4
const REPEATS = 5

L = -0.5 * Matrix{Float64}(I, 8, 8)
problem = IVP(L, ones(8), (0.0, 40.0))
finesolver = RK4(h=1e-5)
coarsesolver = Euler(h=1e-1)

median_elapsed(f; n=REPEATS) = (ts = [@elapsed f() for _ in 1:n]; sort!(ts)[cld(n, 2)])

par = Parareal(finesolver, coarsesolver;
    parameters = PararealParameters(N=N, K=K),
    tolerance  = Tolerance(ϵ=0.0))   # ϵ=0 ⇒ exactly K sweeps, identical work every run

solve(problem, par; mode="SERIAL"); solve(problem, par; mode="THREADS")  # warm up
t_ser = median_elapsed(() -> solve(problem, par; mode="SERIAL"))
t_thr = median_elapsed(() -> solve(problem, par; mode="THREADS"))
P = Threads.nthreads()

# Machine-readable line for the sweep driver to parse, plus a human line.
@printf "RESULT threads=%d N=%d K=%d serial=%.5f threads=%.5f speedup=%.4f\n" P N K t_ser t_thr t_ser/t_thr
