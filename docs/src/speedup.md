# Speedup

Here is a simulation with Parareal on 36 cores, comparing .

```julia
using Distributed, Hwloc
addprocs(31)

@everywhere using RungeKutta, TimeParallel
using BenchmarkTools
using Plots
default(
  framestyle = :box,
  gridalpha = 0.2,
  legend = :none,
  minorgrid = 0.1,
  minorgridstyle = :dash,
  tick_direction = :out
)

u0 = [2.0, 3.0, -14.0]
problem = Lorenz(u0, 0.0, 10.0)

N = 32
tₛ = zeros(N)
tₚ = zeros(N)
for (i, P) in enumerate(1:N)
    solver = Parareal(RK4(h = 1e-4), RK4(h = 1e-2), 𝜑 = 𝜑₂, P = P, K = 1)
    tₛ[i] = @belapsed solve($problem, $solver, mode = "SERIAL")
    tₚ[i] = @belapsed solve($problem, $solver, mode = "DISTRIBUTED")
end

plot(N, tₛ, marker = :o, xscale = :log10)
plot!(N, tₚ, marker = :o, xscale = :log10)
plot!(xlabel = "Number of cores", ylabel = "Time (s)")
```

![svg](imgs/timings.svg)
