@recipe function f(iterate::PararealIterate)
    N = length(iterate)
    for n in 1:N
        @series begin
            if n != 1
                label := ""
            else
                if haskey(plotattributes, :label)
                    label := plotattributes[:label]
                end
            end
            iterate[n]
        end
    end
end

@recipe function f(solution::PararealSolution)
    @↓ lastiterate = solution
    return lastiterate
end

@recipe function f(wrapper::NSDEBase._PhasePlot{<:PararealIterate})
    @↓ iterate ← plottable = wrapper
    N = length(iterate)
    for n in 1:N
        @series begin
            label := if n != 1
                ""
            elseif haskey(plotattributes, :label)
                plotattributes[:label]
            else
                ""
            end
            seriescolor --> 1
            NSDEBase._PhasePlot(iterate[n])
        end
    end
end

@recipe function f(wrapper::NSDEBase._PhasePlot{<:PararealSolution})
    @↓ solution ← plottable = wrapper
    @↓ lastiterate = solution
    return NSDEBase._PhasePlot(lastiterate)
end

@recipe function f(wrapper::NSDEBase._Convergence{<:PararealSolution})
    gridalpha --> 0.2
    markershape --> :circle
    markerstrokewidth --> 0
    seriestype --> :path
    xticks --> 0:1000
    yticks --> 10.0 .^ (-100:100)
    yscale --> :log10
    @↓ solution ← plottable = wrapper
    @↓ errors = solution
    ks = [k for k in eachindex(errors) if isfinite(errors[k]) && errors[k] > 0]
    return ks, errors[ks]
end
