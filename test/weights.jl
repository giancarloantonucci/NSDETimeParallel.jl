@testset "Weight validation" begin
    @test_throws ArgumentError FixedWeight(NaN)
    @test_throws ArgumentError FixedWeight(0)
    for constructor in (SecantWeight, NeighbourWeight)
        @test_throws ArgumentError constructor(w=NaN)
        @test_throws ArgumentError constructor(w=0)
        @test_throws ArgumentError constructor(δ=0)
        @test_throws ArgumentError constructor(δ=Inf)
        @test constructor(w=2).w == 2.0
    end
    @test_throws ArgumentError SecantWeight(τ=-1)
    @test_throws ArgumentError SecantWeight(τ=Inf)
end

@testset "Neighbour measurements, guards, and time units" begin
    U = [[0.0], [1.0], [2.0], [3.0]]
    F = [[0.0], [1.0], [3.0], [7.0]]
    cache = (U=U, F=F, T=[0.0, 1.0, 2.0, 3.0])
    fixed = NSDETimeParallel.WeightState(FixedWeight(2))
    NSDETimeParallel.update_weight!(fixed, cache, nothing, 1)
    @test fixed.w == 2
    state = NSDETimeParallel.WeightState(NeighbourWeight(δ=0.5))
    NSDETimeParallel.update_weight!(state, cache, nothing, 1)
    @test state.w ≈ 8
    NSDETimeParallel.update_weight!(state, cache, nothing, 2)
    @test state.w ≈ 8 # divisor applies to each measurement, not the stored base
    for outputs in (fill([1.0], 4), F)
        flat = merge(cache, (U=fill([1.0], 4), F=outputs))
        NSDETimeParallel.update_weight!(state, flat, nothing, 3)
        @test state.w == 8
    end
    mixed = merge(cache, (U=[[0.0], [1.0], [1.0], [3.0]],))
    state = NSDETimeParallel.WeightState(NeighbourWeight())
    NSDETimeParallel.update_weight!(state, mixed, nothing, 1)
    @test state.w ≈ 2
    doubled = merge(cache, (T=[0.0, 2.0, 4.0, 6.0],))
    state = NSDETimeParallel.WeightState(NeighbourWeight())
    NSDETimeParallel.update_weight!(state, doubled, nothing, 1)
    @test state.w ≈ 2
    λ=0.6
    orbit = (U=[[exp(λ*t)] for t in cache.T], F=[[exp(λ*(t+1))] for t in cache.T], T=cache.T)
    state = NSDETimeParallel.WeightState(NeighbourWeight())
    NSDETimeParallel.update_weight!(state, orbit, nothing, 1)
    @test state.w ≈ exp(2λ) # neighbouring-orbit ratio spans two gains
end

@testset "Same-chunk secant uses latest pair and retains unmeasurable weights" begin
    state = NSDETimeParallel.WeightState(SecantWeight(w=10.0))
    cache = (U=[[1.0], [2.0]], U_=[[0.0], [1.0]], F_=[[0.0], [2.0]], T=[0.0, 1.0, 2.0])
    sol = (lastiterate=[(u=[[2.0]],), (u=[[4.0]],)],)
    NSDETimeParallel.update_weight!(state, cache, sol, 2)
    @test state.w == 2
    @test cache.F_ == [[2.0], [4.0]]
    unchanged = merge(cache, (U_=deepcopy(cache.U),))
    NSDETimeParallel.update_weight!(state, unchanged, sol, 3)
    @test state.w == 2
    @test sol.lastiterate[1].u[end] !== cache.F_[1]
end
