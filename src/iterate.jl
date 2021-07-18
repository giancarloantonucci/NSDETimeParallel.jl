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

getclosest(i₁, t₁, i₂, t₂, target) = (target - t₁ ≥ t₂ - target) ? (i₂, t₂) : (i₁, t₁)

function findclosest(t, target)
    n = length(t)
    if target ≤ t[1]
        return (1, t[1])
    elseif target ≥ t[end]
        return (n, t[end])
    end
    i = 0; j = n; mid = 0
    while i < j
        mid = (i + j) ÷ 2
        if target == t[mid]
            return t[mid]
        elseif target < t[mid]
            if (mid > 1) && (target > t[mid - 1])
                return getclosest(mid - 1, t[mid - 1], mid, t[mid], target)
            end
            j = mid
        else
            if (mid < n) && (target < t[mid + 1])
                return getclosest(mid, t[mid], mid + 1, t[mid + 1], target)
            end
            i = mid + 1
        end
    end
    return (mid, t[mid])
end

function (iterate::TimeParallelIterate)(t::Real)
    N = length(iterate)
    if t ≤ iterate[1].t[1]
        return iterate[1].t[1], iterate[1].u[1]
    elseif t ≥ iterate[end].t[end]
        return iterate[end].t[end], iterate[end].u[end]
    else
        for n = 1:N
            if iterate[n].t[1] ≤ t ≤ iterate[n].t[end]
                (i, tᵢ) = findclosest(iterate[n].t, t)
                uᵢ = iterate[n].u[i]
                return tᵢ, uᵢ
            end
        end
    end
end
