
@testset "Common utilities" begin
  rules_root = joinpath(@__DIR__, "..", "..", "rules")

  # compact rule with weights
  file = joinpath(rules_root, "compact", "CCGV22", "triangle", "4-6.yml")
  qr,data = _read_rule_file(file)

  # expanded rule
  file2 = joinpath(rules_root, "expanded", "WV15", "tetrahedron", "1-1.yml")
  qr2,data2 = _read_rule_file(file2)

  # compact rule without weights
  file3 = joinpath(rules_root, "compact", "Gun75", "7-14.yml")
  qr3,data3 = _read_rule_file(file3)
end