using FiniteElementQuadratureRules
using Documenter
using DocumenterCitations

include(joinpath(@__DIR__, "keylabels.jl"))
include(joinpath(@__DIR__, "generate_rule_overview.jl"))

generate_rule_overview(
    joinpath(@__DIR__, "..", "rules", "compact"),
    joinpath(@__DIR__, "src", "rules"),
)

DocMeta.setdocmeta!(FiniteElementQuadratureRules, :DocTestSetup, :(using FiniteElementQuadratureRules); recursive=true)
bib = CitationBibliography(joinpath(@__DIR__, "..", "references.bib"); style=:keylabels)

makedocs(;
    modules=[FiniteElementQuadratureRules],
    authors="Simon Praetorius <simon.praetorius@tu-dresden.de> and contributors",
    plugins=[bib],
    sitename="FiniteElementQuadratureRules.jl",
    format=Documenter.HTML(;
        canonical="https://spraetor.github.io/FiniteElementQuadratureRules.jl/stable",
        edit_link="main",
        assets=String[],
        size_threshold_warn=200 * 1024,
        size_threshold=512 * 1024,
    ),
    pages=[
        "Home" => "index.md",
        "Database" => [
            "Overview" => "rules/index.md",
            "Point" => "rules/point.md",
            "Line" => "rules/line.md",
            "Triangle" => "rules/triangle.md",
            "Quadrilateral" => "rules/quadrilateral.md",
            "Tetrahedron" => "rules/tetrahedron.md",
            "Hexahedron" => "rules/hexahedron.md",
            "Prism" => "rules/prism.md",
            "Pyramid" => "rules/pyramid.md",
        ],
        "Bibliography" => "bibliography.md",
        "API" => "api.md",
    ],
)

deploydocs(;
    repo="github.com/spraetor/FiniteElementQuadratureRules.jl",
    devbranch="main",
    versions=["stable" => "v^", "v#.#", "dev" => "dev"],
)
