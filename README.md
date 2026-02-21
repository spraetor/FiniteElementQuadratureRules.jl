# FiniteElementQuadratureRules

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://spraetor.github.io/FiniteElementQuadratureRules.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://spraetor.github.io/FiniteElementQuadratureRules.jl/dev/)
[![Build Status](https://github.com/spraetor/FiniteElementQuadratureRules.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/spraetor/FiniteElementQuadratureRules.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/spraetor/FiniteElementQuadratureRules.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/spraetor/FiniteElementQuadratureRules.jl)

This julia package allows to transform compact quadrature rules, given in terms of symmetric orbits of points inside reference domains, like simplices and cubes, into full sets of quadrature points and weights. Additionally, the package provides utilities to improve the accuracy of given quadrature rules, by running an optimizer with high floating point accuracy, starting from the given set of points as initial condition.

## Transforming compact rules
A compat quadrature rule is a set of points given as parameters to symmetric orbits in a reference domain. It can be expanded in a set of points in reference domain coordinates and its associated quadrature weights. The compact format allows to preserve the exact symmetry of the quadrature points, while describing all necessary information to generate the points and weights.

There are some examples of compact rules given in the directory `rules/compact/` extracted from the references listed in the `references.bib` with their corresponding bibtex key.

The first example is a 2d quadrature rule on a triangle domain in `rules/compact/CCGV22/2d/dd2o04_06.yml`

```yml
reference: 'CCGV22'
region: simplex
dim: 2
degree: 4
points: 6
orbits: [0, 2, 0]
arguments:
- '0.10995174365532186763832632490021052896306064753677'
- '0.91576213509770743459571463402201507854325295899829e-1'
- '0.22338158967801146569500700843312280437027268579657'
- '0.44594849091596488631832925388305198839905746639737'
```

This rule can be converted into a full expanded quadrature rules by

```julia
using FiniteElementQuadratureRules
using YAML: load_file, write_file

# read the compact rule from a file and store it in a Dict
data = load_file("rules/compact/CCGV22/2d/dd2o04_06.yml")

# generate the CompactQuadratureRule from the data
cqr = CompactQuadratureRule(data)

# expand the compat rule to generate points and weights
qr = expand(cqr)

# the expanded rule can be written to a .yml file again
mkpath("rules/expanded/CCGV22/2d/")
write_file("rules/expanded/CCGV22/2d/dd2o04_06.yml", Dict(qr, data["reference"]))
```

# Optimizing existing rules
Assume you get a quadrature from an paper, but the printed digits are restricted to low accuracy. For an application you need high accuracy, e.g. usable for quad precision computations. Then you can read an existing rule in compact or expanded form and increase the precision by reducing the quadrature residual on a set of polynomial up to a given degree of exactness. The given points and weights are used as input to the optimizer for an initial condition. Thus, the resulting points should be very close the given input.

We take the rule from above and call an optimizer on its data to improve the accuracy:

```julia
using FiniteElementQuadratureRules
using YAML: load_file

# read the compact rule from a file and store it in a Dict
data = load_file("rules/compact/CCGV22/2d/dd2o04_06.yml")

# generate the CompactQuadratureRule from the data
cqr = CompactQuadratureRule(BigFloat, data)

# now we optimize this rule, resulting in a modified compact rule
oqr = optimize(cqr)
```