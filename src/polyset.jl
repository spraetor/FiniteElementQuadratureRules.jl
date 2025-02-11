# using PolyChaos: AbstractOrthoPoly

struct PolySet{Ω<:AbstractDomain}
  domain::Ω
  basis::Vector{Function}
end

