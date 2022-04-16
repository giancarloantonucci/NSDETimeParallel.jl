"Distributed implementation of Parareal."
function parareal_distributed!(cache::PararealCache, solution::PararealSolution, problem::AbstractInitialValueProblem, parareal::Parareal)
    @↓ iterates, ψ, U, T = solution
    @↓ ℱ, 𝒢, P, K = solver
    @↓ 𝜑, ϵ, Λ, updateΛ = solver.error_check
    # coarse guess
    G = similar(U)
    G[1] = U[1]
    for n = 1:P
        chunk = 𝒢(problem, U[n], T[n], T[n+1])
        G[n+1] = chunk.u[end]
    end
    # main loop
    F = similar(U)
    F[1] = U[1]
    getF(args...) = ℱ(args...).u[end]
    for k = 1:K
        # @↑ solution[k] = U .← U
        solution[k].U .= U
        # for n = 1:k-1
        #     solution[k][n] = solution[k-1][n]
        # end
        # fine run (with Julia's Distributed.jl)
        @sync for n = k:P
            @async F[n+1] = remotecall_fetch(getF, n, problem, U[n], T[n], T[n+1])
        end
        solution[k].F .= F
        # update Lipschitz constant
        Λ = updateΛ ? update_Lipschitz(Λ, U, F) : Λ
        # check convergence
        ψ[k] = 𝜑(solution, k, Λ)
        if ψ[k] ≤ ϵ
            resize!(iterates, k)
            resize!(ψ, k)
            break
        end
        # update (serial)
        for n = k:P
            chunk = 𝒢(problem, U[n], T[n], T[n+1])
            U[n+1] = chunk.u[end] + F[n+1] - G[n+1]
            G[n+1] = chunk.u[end]
        end
        @↑ solution = U, T
    end
    return solution
end
