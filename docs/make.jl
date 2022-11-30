using Documenter
using TimeParallel

PAGES = ["Home" => "index.md"]

makedocs(;
    sitename = "TimeParallel.jl",
    format = Documenter.HTML(),
    modules = [TimeParallel],
    pages = PAGES,
    authors = "Giancarlo A. Antonucci <giancarlo.antonucci@icloud.com>"
)

deploydocs(;
    repo = "https://github.com/giancarloantonucci/TimeParallel.jl"
)
