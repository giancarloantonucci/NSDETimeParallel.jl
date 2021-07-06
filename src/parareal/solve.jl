function solve!(solution::TimeParallelSolution, problem, solver::Parareal)
    @↓ iterates, φ, U, T = solution
    @↓ rhs = problem
    @↓ mode, ℱ, 𝒢, P, K = solver
    @↓ 𝜑, ϵ = solver.objective

    G = similar(U)
    G[1] = U[1]
    for n = 1:P
        chunk = 𝒢(rhs, U[n], T[n], T[n+1])
        G[n+1] = chunk.u[end]
    end
    # G .= U

    k = 0
    F = similar(U); F[1] = U[1]
    for outer k = 1:K
        # fine run (parallelisable)
        # if mode == "SERIAL"
        @↑ solution[k] = U, T
        for n = k:P
            chunk = ℱ(rhs, U[n], T[n], T[n+1])
            solution[k][n] = chunk
            F[n+1] = chunk.u[end]
        end
        # fine run (uses Julia's Distributed)
        # elseif mode == "DISTRIBUTED"
        #     v = pmap(n -> ℱ(rhs, U[n], T[n], T[n+1]), 1:P)
        #     for n = k:P
        #         chunk = v[n]
        #         solution[k][n] = chunk
        #         F[n+1] = chunk.u[end]
        #     end
        # end
        # check convergence
        φ[k] = 𝜑(U, F, T)
        if φ[k] ≤ ϵ
            break
        end
        # update (serial)
        for n = k:P
            chunk = 𝒢(rhs, U[n], T[n], T[n+1])
            U[n+1] = chunk.u[end] + F[n+1] - G[n+1]
            G[n+1] = chunk.u[end]
        end
        @↑ solution = U, T
    end
    resize!(iterates, k)
    resize!(φ, k)
    return solution
end
