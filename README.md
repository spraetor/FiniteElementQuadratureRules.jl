# FiniteElementQuadratureRules

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://spraetor.github.io/FiniteElementQuadratureRules.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://spraetor.github.io/FiniteElementQuadratureRules.jl/dev/)
[![Build Status](https://github.com/spraetor/FiniteElementQuadratureRules.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/spraetor/FiniteElementQuadratureRules.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/spraetor/FiniteElementQuadratureRules.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/spraetor/FiniteElementQuadratureRules.jl)

This julia package allows to transform compact quadrature rules, given in terms of symmetric orbits of points inside reference domains, like simplices and cubes, into full sets of quadrature points and weights. Additionally, the package provides utilities to improve the accuracy of given quadrature rules, by running an optimizer with high floating point accuracy, starting from the given set of points as initial condition.

## Main workflow for FE library developers
The main purpose of this repository is to help finite-element library developers maintain quadrature-rule tables in their own code base:

1. Read compact rules from `rules/compact/...`.
2. Expand them to full coordinates and weights.
3. Transform those rules to the reference-element convention used by your FE library.
4. Generate library-specific source/header files from a template.

The package includes a generator template for Dune: `dune.templ.hh`.

## Expand compact rules
A compact rule encodes symmetric orbits plus orbit parameters. Expanding it yields explicit quadrature coordinates and weights.

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
Different FE libraries use different reference-element coordinate systems. Use `transform` to map one convention to another.

Example (transforming to Dune coordinates):

```julia
using FiniteElementQuadratureRules
using YAML: load_file

data = load_file("rules/compact/CCGV22/triangle/4-6.yml")
qr = expand(CompactQuadratureRule(Float64, data))

ref_dune = duneReferenceElement(domain(qr))
qr_dune = transform(qr, ref_dune)
```

## Generate FE-library specific code from templates
The `generate` function reads compact-rule files, expands/sorts/selects rules, and renders output files from an Otera template.

Example for Dune headers:

```julia
using FiniteElementQuadratureRules

generate("dune.templ.hh", "rules/compact/CCGV22/", "dune/"; precision=80)
```

This creates one output file per domain in `dune/` (e.g. `triangle.hh`, `tetrahedron.hh`).

## Optional: improve rule accuracy
Accuracy optimization is an additional feature. It is useful when a published rule has limited printed precision and you need a high-precision variant.

```julia
using FiniteElementQuadratureRules
using YAML: load_file

data = load_file("rules/compact/CCGV22/triangle/4-6.yml")
cqr = CompactQuadratureRule(BigFloat, data)
oqr = optimize(cqr)
```
