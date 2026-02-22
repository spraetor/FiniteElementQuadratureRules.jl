using FiniteElementQuadratureRules
using Documenter

DocMeta.setdocmeta!(FiniteElementQuadratureRules, :DocTestSetup, :(using FiniteElementQuadratureRules); recursive=true)

makedocs(;
    modules=[FiniteElementQuadratureRules],
    authors="Simon Praetorius <simon.praetorius@tu-dresden.de> and contributors",
    sitename="FiniteElementQuadratureRules.jl",
    format=Documenter.HTML(;
        canonical="https://spraetor.github.io/FiniteElementQuadratureRules.jl/stable",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/spraetor/FiniteElementQuadratureRules.jl",
    devbranch="main",
    versions=["stable" => "v^", "v#.#", "dev" => "dev"],
)
