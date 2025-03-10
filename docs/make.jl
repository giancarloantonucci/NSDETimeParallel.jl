using Documenter
using DocumenterInterLinks
using NSDEBase, NSDETimeParallel

PAGES = ["Home" => "index.md"]

links = InterLinks(
    "NSDEBase" => (
        "https://giancarloantonucci.github.io/NSDEBase.jl/dev/",
        "https://giancarloantonucci.github.io/NSDEBase.jl/dev/objects.inv"
    )
)

makedocs(;
    sitename = "NSDETimeParallel.jl",
    format = Documenter.HTML(),
    modules = [NSDETimeParallel],
    pages = PAGES,
    authors = "Giancarlo A. Antonucci <giancarlo.antonucci@icloud.com>",
    plugins = [links],
)

deploydocs(;
    repo = "https://github.com/giancarloantonucci/NSDETimeParallel.jl"
)
