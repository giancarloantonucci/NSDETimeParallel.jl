using Documenter
using NSDETimeParallel

PAGES = [
    "Home" => "index.md",
    "Conventions" => "conventions.md",
    "Backends" => "backends.md",
    "Criterion" => "criterion.md",
    "API" => "api.md"
]

makedocs(;
    sitename = "NSDETimeParallel.jl",
    format = Documenter.HTML(),
    modules = [NSDETimeParallel],
    pages = PAGES,
    checkdocs = :exports, # every export must carry a docstring, or the build fails
    authors = "Giancarlo A. Antonucci <giancarlo.antonucci@icloud.com>"
)

deploydocs(;
    repo = "github.com/giancarloantonucci/NSDETimeParallel.jl.git"
)
