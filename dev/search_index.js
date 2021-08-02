var documenterSearchIndex = {"docs":
[{"location":"#TimeParallel.jl","page":"Home","title":"TimeParallel.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"#Public-API","page":"Home","title":"Public API","text":"","category":"section"},{"location":"#Constructors","page":"Home","title":"Constructors","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"solve\nsolve!\nTimeParallelSolver\nTimeParallelIterate\nTimeParallelSolution","category":"page"},{"location":"#NSDEBase.solve","page":"Home","title":"NSDEBase.solve","text":"solve(problem, solver::TimeParallelSolver; mode::String = \"SERIAL\")\n\nreturns the TimeParallelSolution of a problem, e.g. an InitialValueProblem.\n\n\n\n\n\n","category":"function"},{"location":"#TimeParallel.TimeParallelSolution","page":"Home","title":"TimeParallel.TimeParallelSolution","text":"TimeParallelSolution{iterates_T, φ_T, U_T, T_T} <: InitialValueSolution\n\nreturns a constructor for the numerical solution of an InitialValueProblem obtained with a TimeParallelSolver.\n\n\n\nTimeParallelSolution(problem, solver::TimeParallelSolver)\n\nreturns an initialised TimeParallelSolution given a problem, e.g. an InitialValueProblem, and a TimeParallelSolver.\n\n\n\n\n\n","category":"type"},{"location":"#Solvers","page":"Home","title":"Solvers","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Parareal","category":"page"},{"location":"#TimeParallel.Parareal","page":"Home","title":"TimeParallel.Parareal","text":"Parareal{error_check_T, ℱ_T, 𝒢_T, P_T, K_T} <: TimeParallelSolver\n\nreturns a constructor for the TimeParallelSolver based on the parareal algorithm.\n\n\n\nParareal(ℱ::Function, 𝒢::Function; P = 10, K = P, 𝜑 = 𝜑₁, ϵ = 1e-12)\n\nreturns a Parareal with:\n\nℱ :: Function : fine solver.\n𝒢 :: Function : coarse solver.\nP :: Integer  : number of time chunks.\nK :: Integer  : maximum number of iterations.\n𝜑 :: Function : error control function.\nϵ :: Real     : tolerance.\n\n\n\nParareal(finesolver::InitialValueSolver, coarsolver::InitialValueSolver; P = 10, K = P, 𝜑 = 𝜑₁, ϵ = 1e-12)\n\nreturns a Parareal from a finesolver and a coarsolver.\n\n\n\n\n\n","category":"type"},{"location":"#Utilities","page":"Home","title":"Utilities","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"ErrorCheck\ncoarseguess!\ngetchunks","category":"page"},{"location":"#Index","page":"Home","title":"Index","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"}]
}
