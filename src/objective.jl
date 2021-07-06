struct ErrorFunction{𝜑_T, ϵ_T}
    𝜑::𝜑_T
    ϵ::ϵ_T
end

function 𝜑₁(U, F, T)
    r = 0.0
    for n = 2:length(U)
        r += norm(U[n] - F[n]) / norm(U[n])
    end
    return r
end

function 𝜑₂(U, F, T, λ = 0.0)
    r = 0.0
    for n = 2:length(U)
        Wₙ = exp(λ * (T[1] - T[n]))
        r += norm(Wₙ * (U[n] - F[n]))^2
    end
    return r
end
