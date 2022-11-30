using NSDETimeParallel
using Documenter

DocMeta.setdocmeta!(NSDETimeParallel, :DocTestSetup, :(using NSDETimeParallel); recursive=true)

makedocs(;
    modules=[NSDETimeParallel],
    authors="Giancarlo A. Antonucci",
    repo="https://github.com/giancarloantonucci/NSDETimeParallel.jl/blob/{commit}{path}#{line}",
    sitename="NSDETimeParallel.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://giancarloantonucci.github.io/NSDETimeParallel.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/giancarloantonucci/NSDETimeParallel.jl",
    devbranch="main",
)
