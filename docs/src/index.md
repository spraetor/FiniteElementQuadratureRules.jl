```@meta
CurrentModule = FiniteElementQuadratureRules
```

# FiniteElementQuadratureRules

This julia package allows to transform compact quadrature rules, given in terms of symmetric orbits of points inside reference domains, like simplices and cubes, into full sets of quadrature points and weights. Additionally, the package provides utilities to improve the accuracy of given quadrature rules, by running an optimizer with high floating point accuracy, starting from the given set of points as initial condition.

## Main workflow for FE library developers
The primary goal is to support finite-element library developers who maintain stored quadrature rules:

1. read compact rules from `rules/compact/...`,
2. expand to full coordinates and weights,
3. transform coordinates/weights to the target reference-element convention,
4. generate library-specific code from a template.

The repository includes a template for Dune in `dune.templ.hh`.

## Expand compact rules
A compact rule describes symmetric orbits and orbit parameters. Expansion yields explicit quadrature coordinates and weights.

Example:

```julia
using FiniteElementQuadratureRules
using YAML: load_file, write_file

data = load_file("rules/compact/CCGV22/triangle/4-6.yml")
cqr = CompactQuadratureRule(BigFloat, data)
qr = expand(cqr)

mkpath("rules/expanded/CCGV22/triangle")
write_file("rules/expanded/CCGV22/triangle/4-6.yml", Dict(qr; reference=data["reference"]))
```

## Transform to another reference-element convention
Finite-element libraries may define different reference-element coordinates for the same topological domain. Use `transform` to map rules between conventions.

Example (to Dune triangle convention):

```julia
using FiniteElementQuadratureRules
using StaticArrays: SVector
using YAML: load_file

data = load_file("rules/compact/CCGV22/triangle/4-6.yml")
qr = expand(CompactQuadratureRule(Float64, data))

ref_in_int = ReferenceElement(Triangle())
ref_in = ReferenceElement{Triangle,SVector{2,Float64}}(
  [SVector{2,Float64}(x) for x in ref_in_int.coordinates],
  ref_in_int.facets
)
ref_dune_int = duneReferenceElement(Triangle())
ref_dune = ReferenceElement{Triangle,SVector{2,Float64}}(
  [SVector{2,Float64}(x) for x in ref_dune_int.coordinates],
  ref_dune_int.facets
)

qr_dune = transform(qr, ref_in, ref_dune)
```

## Generate library-specific code from templates
`generate` renders source/header files from an Otera template with rule data grouped by domain and degree.

Example (Dune):

```julia
using FiniteElementQuadratureRules

generate("dune.templ.hh", "rules/compact/CCGV22/", "dune/"; precision=80)
```

## Optional: optimize rule accuracy
Improving rule accuracy is an additional feature. Use it when published coefficients are low precision and you need a high-precision variant.

```julia
using FiniteElementQuadratureRules
using YAML: load_file

data = load_file("rules/compact/CCGV22/triangle/4-6.yml")
cqr = CompactQuadratureRule(BigFloat, data)
oqr = optimize(cqr)
```
