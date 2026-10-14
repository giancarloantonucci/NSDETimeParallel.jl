"A stopping rule evaluated on the fine trajectories from the current starts."
abstract type AbstractCriterion end
"All-interface mixed update test, with explicit absolute and relative tolerances."
struct RelativeUpdate{T<:Real} <: AbstractCriterion
    atol::T
    rtol::T
    function RelativeUpdate(; atol::Real=1e-12, rtol::Real=1e-9)
        isfinite(atol) && isfinite(rtol) && atol >= 0 && rtol >= 0 && atol + rtol > 0 ||
            throw(ArgumentError("finite nonnegative tolerances, at least one positive, required"))
        a, r = promote(float(atol), float(rtol))
        return new{typeof(a)}(a, r)
    end
end
"Normalized weighted defect p-mean; p=1 is the average and p=Inf the maximum."
struct WeightedDefect{T<:Real,R<:AbstractWeightRule,P<:Real} <: AbstractCriterion
    ϵ::T
    weight::R
    p::P
end
function WeightedDefect(; ϵ::Real=1e-9, weight=FixedWeight(), p::Real=1)
    isfinite(ϵ) && ϵ >= 0 || throw(ArgumentError("tolerance must be finite and nonnegative"))
    p >= 1 || throw(ArgumentError("p must be at least one"))
    return WeightedDefect(float(ϵ), weight_rule(weight), p)
end
"Legacy mean relative interface-update criterion."
struct MeanRelativeUpdate{T<:Real} <: AbstractCriterion
    ϵ::T
    function MeanRelativeUpdate(; ϵ::Real=1e-12)
        isfinite(ϵ) && ϵ >= 0 || throw(ArgumentError("tolerance must be finite and nonnegative"))
        return new{typeof(float(ϵ))}(float(ϵ))
    end
end
"Disable tolerance acceptance while retaining budgets, failure checks, and exact completion."
struct NoStop <: AbstractCriterion end
weight_rule(::AbstractCriterion) = FixedWeight()
weight_rule(c::WeightedDefect) = c.weight
threshold(c::AbstractCriterion) = c.ϵ
threshold(::RelativeUpdate) = 1.0
threshold(::NoStop) = 0.0
"Return whether a criterion accepts its evaluated residual."
accept(c::AbstractCriterion, value) = isfinite(value) && value <= threshold(c)
accept(::NoStop, value) = false
"Evaluate a criterion from the current fine batch and its cached starts."
residual(::NoStop, cache, k) = Inf
function residual(c::RelativeUpdate, cache, k)
    k == 1 && return Inf
    value = zero(c.atol)
    for n in eachindex(cache.U)
        u, v = cache.U[n], cache.U_[n]
        change = norm(u - v)
        scale = c.atol + c.rtol * max(norm(u), norm(v))
        q = iszero(scale) ? (iszero(change) ? zero(value) : Inf) : change / scale
        isfinite(q) || return Inf
        value = max(value, q)
    end
    return value
end
function residual(::MeanRelativeUpdate, cache, k)
    previous = k > 1 ? cache.U_ : cache.F
    value = 0.0
    for n in eachindex(cache.U)
        scale = norm(cache.U[n])
        change = norm(cache.U[n]-previous[n])
        value += iszero(scale) ? change : change/scale
    end
    return value / length(cache.U)
end
function residual(c::WeightedDefect, cache, k)
    w = max(one(cache.weightstate.w), cache.weightstate.w)
    value = 0.0
    if c.p == 1 || c.p == Inf
        for n in eachindex(cache.U)
            defect = norm(cache.U[n]-cache.F[n])*w^(cache.T[1]-cache.T[n])
            value = c.p == 1 ? value+defect : max(value, defect)
        end
        return c.p == 1 ? value/length(cache.U) : value
    end
    values = [norm(cache.U[n]-cache.F[n])*w^(cache.T[1]-cache.T[n]) for n in eachindex(cache.U)]
    return norm(values, c.p)/length(values)^(1/c.p)
end

"Pipeline-only relative frontier test. It is not the bulk all-interface update test."
struct FrontierRelativeUpdate{T<:Real} <: AbstractCriterion
    ϵ::T
    function FrontierRelativeUpdate(; ϵ::Real=1e-9)
        isfinite(ϵ) && ϵ >= 0 || throw(ArgumentError("tolerance must be finite and nonnegative"))
        return new{typeof(float(ϵ))}(float(ϵ))
    end
end
function residual(::FrontierRelativeUpdate, cache, k)
    return throw(ArgumentError("FrontierRelativeUpdate requires PipelinedMPIBackend"))
end
function pipeline_kind(::AbstractCriterion)
    return throw(
        ArgumentError(
            "criterion is not supported by PipelinedMPIBackend; use a fixed-weight maximum defect, FrontierRelativeUpdate, or NoStop",
        ),
    )
end
pipeline_kind(::NoStop) = :none
pipeline_kind(::FrontierRelativeUpdate) = :frontier
function pipeline_kind(c::WeightedDefect)
    c.p == Inf && c.weight isa FixedWeight ||
        throw(ArgumentError("pipeline weighted acceptance requires p=Inf and a fixed weight"))
    return :defect
end
