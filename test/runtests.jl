using NSDETimeParallel
using Test

using Pkg
for (uuid, pkg) in Pkg.dependencies()
    if pkg.name == "NSDEBase"
        @info "Testing NSDERungeKutta with NSDEBase version: $(pkg.version) from $(pkg.source)"
    end
end

@testset "ErrorControl" begin
    # Write your tests here.
end

@testset "Coarse" begin
    # Write your tests here.
end

@testset "Parareal" begin
    # Write your tests here.
end
