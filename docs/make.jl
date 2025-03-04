using UrbanTethysChloris
using Documenter

DocMeta.setdocmeta!(
    UrbanTethysChloris,
    :DocTestSetup,
    :(using UrbanTethysChloris);
    recursive = true,
)

const page_rename = Dict("developer.md" => "Developer docs") # Without the numbers
const numbered_pages = [
    file for file in readdir(joinpath(@__DIR__, "src")) if
    file != "index.md" && splitext(file)[2] == ".md"
]

makedocs(;
    modules = [UrbanTethysChloris],
    authors = "Hugo Solleder <hugo.solleder@epfl.ch>",
    repo = "https://github.com/EPFL-ENAC/UrbanTethysChloris.jl/blob/{commit}{path}#{line}",
    sitename = "UrbanTethysChloris.jl",
    format = Documenter.HTML(;
        canonical = "https://EPFL-ENAC.github.io/UrbanTethysChloris.jl",
    ),
    pages = ["index.md"; numbered_pages],
)

deploydocs(; repo = "github.com/EPFL-ENAC/UrbanTethysChloris.jl")
