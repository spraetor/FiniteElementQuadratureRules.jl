module FiniteElementQuadratureRulesExportExt

using Base: Filesystem
using Printf: @sprintf
using FiniteElementQuadratureRules
using OteraEngine
import BibParser
import BibInternal
import BibFormatter
import YAML

import FiniteElementQuadratureRules: default_chooser, duneReferenceElement, generate, generate_rule_overview

include("dune.jl")
include("generate_common.jl")
include("generate.jl")
include("generate_rule_overview.jl")

end
