"An abstract type for time-parallel solvers."
abstract type TimeParallelSolver <: InitialValueSolver end

# ---------------------------------------------------------------------------- #
#                                   Functions                                  #
# ---------------------------------------------------------------------------- #

"""
    show(io::IO, solver::TimeParallelSolver)

prints a full description of `solver` and its contents to a stream `io`.
"""
Base.show(io::IO, solver::TimeParallelSolver) = NSDEBase._show(io, solver)

"""
    summary(io::IO, solver::TimeParallelSolver)

prints a brief description of `solver` to a stream `io`.
"""
Base.summary(io::IO, solver::TimeParallelSolver) = NSDEBase._summary(io, solver)
