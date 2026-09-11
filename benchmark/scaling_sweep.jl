# Parareal pool-scaling curve: launches benchmark/scaling.jl at a range of
# thread counts (the pool cannot be resized within a live Julia session) and
# collates the speed-up against pool size.
#
#   julia --project=benchmark benchmark/scaling_sweep.jl [N] [K] [maxpool]
#
# The pool is swept 1, 2, 4, 8, … up to maxpool (default: this machine's core
# count). Speed-up should climb until pool ≈ N, then flatten as extra workers
# find no chunk to take, and is bounded above by the serial coarse+correction.
using Printf

N       = length(ARGS) ≥ 1 ? parse(Int, ARGS[1]) : 16
K       = length(ARGS) ≥ 2 ? parse(Int, ARGS[2]) : 4
maxpool = length(ARGS) ≥ 3 ? parse(Int, ARGS[3]) : Sys.CPU_THREADS

script  = joinpath(@__DIR__, "scaling.jl")
project = dirname(Base.active_project()) # the DIR, not the Project.toml path — child julia needs the folder
jl      = Base.julia_cmd()[1]

pools = [2^k for k in 0:floor(Int, log2(maxpool))]
pools[end] == maxpool || push!(pools, maxpool)

@printf "Parareal pool scaling — N=%d, K=%d (fixed work), fine RK4(h=1e-5), 8-dim linear decay.\n" N K
@printf "Sweeping pool over %s on a %d-core machine.\n\n" pools Sys.CPU_THREADS
@printf "%6s %11s %12s %10s %12s\n" "pool" "serial(s)" "threads(s)" "speedup" "efficiency"

base_speedup = nothing
for pool in pools
    out = read(`$jl --project=$project --threads=$pool $script $N $K`, String)
    m = match(r"serial=([\d.]+) threads=([\d.]+) speedup=([\d.]+)", out)
    m === nothing && (println("  (pool=$pool: no result — ", strip(out), ")"); continue)
    ser, thr, sp = parse.(Float64, m.captures)
    eff = sp / pool                     # speed-up per worker: 1.0 = perfect
    @printf "%6d %11.3f %12.3f %10.2f %11.0f%%\n" pool ser thr sp 100eff
end

println("\nSpeed-up climbs with the pool until pool ≈ N, then flattens (idle workers).")
println("The plateau height is the Amdahl ceiling: only the fine sweep is parallel.")
