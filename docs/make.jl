using Documenter
using TimeParallel

PAGES = ["Home" => "index.md"]

makedocs(;
    sitename = "TimeParallel",
    format = Documenter.HTML(),
    modules = [TimeParallel],
    pages = PAGES,
    authors = "Giancarlo A. Antonucci <giancarlo.antonucci@icloud.com>"
)

deploydocs(;
    repo = "https://github.com/antonuccig/TimeParallel.jl"
)
