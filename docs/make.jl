using Documenter
using NSDETimeParallel

PAGES = [
    "Home" => "index.md",
    "Guided tour" => "guide.md",
    "Conventions" => "conventions.md",
    "Backends" => "backends.md",
    "Criterion" => "criterion.md",
    "Extensions" => "extensions.md",
    "Quick reference" => "quick_reference.md",
    "API" => "api.md",
]

makedocs(;
    root=@__DIR__,
    sitename="NSDETimeParallel.jl",
    format=Documenter.HTML(),
    modules=[NSDETimeParallel],
    pages=PAGES,
    checkdocs=:exports, # every export must carry a docstring, or the build fails
    authors="Giancarlo A. Antonucci <giancarlo.antonucci@icloud.com>",
)

deploydocs(; repo="github.com/giancarloantonucci/NSDETimeParallel.jl.git")
