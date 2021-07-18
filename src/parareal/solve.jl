function solve!(solution::TimeParallelSolution, problem, solver::Parareal)
    @↓ iterates, φ, U, T = solution
    @↓ mode, ℱ, 𝒢, P, K = solver
    @↓ 𝜑, ϵ = solver.objective
    # coarse guess
    G = similar(U); G[1] = U[1]
    for n = 1:P
        chunk = 𝒢(problem, U[n], T[n], T[n+1])
        G[n+1] = chunk.u[end]
    end
    # main loop
    F = similar(U); F[1] = U[1]
    @everywhere problem = $problem
    for k = 1:K
        # fine run (parallelisable)
        if nprocs() == 1 || mode == "SERIAL"
            for n = k:P
                chunk = ℱ(problem, U[n], T[n], T[n+1])
                solution[k][n] = chunk
                F[n+1] = chunk.u[end]
            end
        elseif nprocs() > 1 && P == nworkers()
            ws = (k:P)
            @everywhere ws begin
                n = myid()
                U = $U
                T = $T
                Un = U[n]
                Tn = T[n]
                Tp = T[n+1]
                chunk = ℱ(problem, Un, Tn, Tp)
            end
            for n in ws
                chunk = @fetchfrom n Main.chunk
                solution[k][n] = chunk
                F[n+1] = chunk.u[end]
            end
        end
        # check convergence
        φ[k] = 𝜑(U, F, T)
        if φ[k] ≤ ϵ
            resize!(iterates, k)
            resize!(φ, k)
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
