mutable struct TimeParallelIterate{chunks_T, U_T, T_T}
    chunks::chunks_T
    U::U_T
    T::T_T
end

function TimeParallelIterate(problem, solver::TimeParallelSolver)
    @↓ u0, t0 ← tspan[1] = problem
    @↓ P = solver
    chunks = Vector{Any}(undef, P)
    U = Vector{typeof(u0)}(undef, P+1)
    T = Vector{typeof(t0)}(undef, P+1)
    TimeParallelIterate(chunks, U, T)
end

# solution[k][n] ≡ solution.iterates[k].chunks[n]
Base.length(iterate::TimeParallelIterate) = length(iterate.chunks)
Base.getindex(iterate::TimeParallelIterate, n::Int) = iterate.chunks[n]
Base.setindex!(iterate::TimeParallelIterate, value, n::Int) = iterate.chunks[n] = value
Base.lastindex(iterate::TimeParallelIterate) = lastindex(iterate.chunks)

# ----------------------------------- MISC ----------------------------------- #

function (iterate::TimeParallelIterate)(t::Real)
    N = length(iterate)
    if t ≤ iterate[1].t[1]
        return iterate[1].t[1], iterate[1].u[1]
    elseif t ≥ iterate[end].t[end]
        return iterate[end].t[end], iterate[end].u[end]
    else
        for n = 1:N
            if iterate[n].t[1] ≤ t ≤ iterate[n].t[end]
                return iterate[n](t)
            end
        end
    end
end
