using Printf

function default_chooser(qr1_entry, qr2_entry)
  qr1, qr2 = qr1_entry[1], qr2_entry[1]
  if qr1.degree != qr2.degree
    return qr1.degree < qr2.degree
  elseif length(qr1) != length(qr2)
    return length(qr1) < length(qr2)
  elseif (:inside in qr1.properties) != (:inside in qr2.properties)
    return :inside in qr1.properties
  elseif (:positive in qr1.properties) != (:positive in qr2.properties)
    return :positive in qr1.properties
  else
    entry1, entry2 = qr1_entry[2], qr2_entry[2]
    if length(entry2.date.year) == 0
      return true
    elseif length(entry1.date.year) == 0
      return false
    else
      year1 = parse(Int, entry1.date.year)
      year2 = parse(Int, entry2.date.year)
      return year1 < year2
    end
  end
end


function _read_rule_file(path::AbstractString)
  data = YAML.load_file(path)
  qr = if haskey(data, "coordinates")
    QuadratureRule(BigFloat, data)
  elseif haskey(data, "weights")
    expand(CompactQuadratureRuleWithWeights(BigFloat, data))
  else
    expand(CompactQuadratureRule(BigFloat, data))
  end
  return data, qr
end


_markdown_escape(s::AbstractString) = replace(s, "|" => "\\|", "\n" => " ")
_markdown_escape(x) = _markdown_escape(string(x))

function _format_accuracy(x)
  value = x isa AbstractString ? parse(BigFloat, x) : BigFloat(x)
  return @sprintf("%.2e", Float64(value))
end


function _format_rule_reference(reference::AbstractString)
  return "[@$(_markdown_escape(reference))]"
end


function _format_orbits(data::Dict)
  haskey(data, "orbits") || return "-"
  "[" * join(string.(data["orbits"]), ", ") * "]"
end
