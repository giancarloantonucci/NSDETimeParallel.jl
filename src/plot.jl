@recipe function f(iterate::TimeParallelIterate; vars = nothing, chunks = false, label = "")
    fontfamily --> "Computer Modern"
    framestyle --> :box
    gridalpha --> 0.2
    linewidth --> 1.5
    minorgrid --> 0.1
    minorgridstyle --> :dash
    seriestype --> :path
    tick_direction --> :out
    xwiden --> false
    P = length(iterate)
    for n in 1:P
        @series begin
            if n != P
                label := ""
            else
                label := label
            end
            vars --> vars
            seriescolor --> n
            iterate[n]
        end
    end
    if chunks == true
        @series begin
            label := ""
            seriescolor := :black
            seriestype := :vline
            linealpha --> 0.5
            linestyle --> :dashdotdot
            getchunks(iterate)
        end
    end
    # solves a bug in Plots
    primary := false
    ()
end

@recipe function f(solution::TimeParallelSolution; vars = nothing, chunks = false, label = "")
    @series begin
        vars   --> vars
        chunks --> chunks
        solution[end]
    end
end

@userplot CONVERGENCE
@recipe function f(h::CONVERGENCE)
    if length(h.args) != 1 || !(h.args[1] isa TimeParallelSolution)
        error("convergence must be given TimeParallelSolution. Got $(typeof(h.args))")
    elseif h.args[1] isa TimeParallelSolution
        solution = h.args[1]
    end
    fontfamily      --> "Computer Modern"
    framestyle      --> :box
    legend          --> :none
    gridalpha       --> 0.2
    linewidth       --> 1.5
    markershape     --> :circle
    markersize      --> 3
    seriestype      --> :path
    tick_direction  --> :out
    yminorgrid      --> true
    yminorgridstyle --> :dash
    solution.φ
end
