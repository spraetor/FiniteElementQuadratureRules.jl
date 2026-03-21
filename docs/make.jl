using FiniteElementQuadratureRules
using Documenter
using DocumenterCitations

DocMeta.setdocmeta!(FiniteElementQuadratureRules, :DocTestSetup, :(using FiniteElementQuadratureRules); recursive=true)
bib = CitationBibliography(joinpath(@__DIR__, "..", "references.bib"))

makedocs(;
    modules=[FiniteElementQuadratureRules],
    authors="Simon Praetorius <simon.praetorius@tu-dresden.de> and contributors",
    plugins=[bib],
    sitename="FiniteElementQuadratureRules.jl",
    format=Documenter.HTML(;
        canonical="https://spraetor.github.io/FiniteElementQuadratureRules.jl/stable",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "All Rules" => "allrules.md",
        "API" => "api.md",
    ],
)

deploydocs(;
    repo="github.com/spraetor/FiniteElementQuadratureRules.jl",
    devbranch="main",
    versions=["stable" => "v^", "v#.#", "dev" => "dev"],
)
