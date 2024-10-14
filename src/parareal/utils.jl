function theoretical_speedup(K::Integer, N::Integer, ζ::AbstractFloat)
    term1 = 1 + ζ * N
    term2 = 1 - (K-1) / 2N
    return N / (K * term1 * term2)
end
