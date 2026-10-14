@testset "Bulk accounting and returned state" begin
    p=Logistic([0.1], (0.0, 1.0))
    F=RK4(h=0.01)
    G=Euler(h=0.05)
    for N in (1, 4), mode in (SerialBackend(), ThreadsBackend())
        a=Parareal(F, G; parameters=PararealParameters(N=N), criterion=NoStop())
        s=solve(p, a; backend=mode, history=:interfaces)
        d=s.diagnostics
        @test d.status==:exact
        @test d.fine_batches==N
        @test d.sweeps==N-1
        @test d.fine_chunk_calls==N*(N+1)÷2
        @test d.coarse_chunk_calls==N*(N-1)÷2
        u=boundaryvalues(s.lastiterate)
        truth=solve(p, F)
        @test u[end]==truth.u[end]
        @test maxseam(s.lastiterate)==0
        @test length(d.interfaces)==N
    end
    a=Parareal(F, G; parameters=PararealParameters(N=4, K=2), criterion=NoStop())
    c=initialize_cache(p, a)
    s=initialize_solution(p, a)
    a(c, s, p; history=:interfaces)
    @test s.diagnostics.status==:budget_exhausted
    @test s.diagnostics.sweeps==1
    @test boundaryvalues(s.lastiterate)==s.diagnostics.interfaces[end]
    first=deepcopy(boundaryvalues(s.lastiterate))
    a(c, s, p; history=:interfaces)
    @test boundaryvalues(s.lastiterate)==first
end

@testset "nonfinite coarse pass has no fabricated fine output" begin
    p = Dahlquist([1.0], (0.0, 8.0); λ=-1000.0)
    bad = Parareal(RK4(h=0.001), Euler(h=0.02); N=8)
    for backend in (SerialBackend(), ThreadsBackend())
        cache = initialize_cache(p, bad)
        buffer = initialize_solution(p, bad)
        for _ in 1:2
            result = bad(cache, buffer, p; backend)
            @test result.diagnostics.status == :nonfinite
            @test result.diagnostics.fine_batches == 0
            @test isempty(result.errors)
            @test all(u -> all(isnan, u), boundaryvalues(result.lastiterate))
            @test all(t -> all(isnan, result(t)), (0.0, 0.5, 4.0, 8.0))
            @test boundarytimes(result.lastiterate) == collect(0.0:8.0)
            @test_throws ArgumentError Wnorm(result, reference(p, bad))
            # Reuse buffers filled with NaN, then fail again after real output.
            good = remake(bad; coarsesolver=RK4(h=0.001))
            good(cache, buffer, p; backend)
            @test all(u -> all(isfinite, u), boundaryvalues(buffer.lastiterate))
        end
    end
end
