@recipe function f(iterate::PararealIterate; label="")
    P = length(iterate)
    for n = 1:P
        @series begin
            label := n != P ? "" : label
            iterate[n]
        end
    end
    # solves a bug in Plots with OrdinaryDiffEq
    # primary := false
    # ()
end

@recipe function f(solution::PararealSolution)
    @↓ saveiterates = solution
    @series begin
        if saveiterates
            @↓ alliterates = solution
            return alliterates[end]
        else
            @↓ lastiterate = solution
            return lastiterate
        end
    end
end

@userplot CONVERGENCE
@recipe function f(h::CONVERGENCE; label="")
    if length(h.args) == 1 && h.args[1] isa PararealSolution
        solution = h.args[1]
    else
        error("convergence must be given PararealSolution. Got $(typeof(h.args))")
    end
    framestyle --> :box
    gridalpha --> 0.2
    linewidth --> 1.5
    markershape --> :circle
    markersize --> 3
    minorgrid --> 0.1
    minorgridstyle --> :dash
    seriestype --> :path
    tick_direction --> :out
    @series begin
        label := label
        solution.errors
    end
end
