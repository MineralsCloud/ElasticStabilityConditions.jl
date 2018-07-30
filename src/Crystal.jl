"""
# module Crystal

- Julia version: 0.7
- Author: qz
- Date: 2018-07-25

# Examples

```jldoctest
julia>
```
"""
module Crystal

struct CrystalSystem
    elastic_matrix::Matrix
    CrystalSystem(elastic_matrix) = size(elastic_matrix) == (6, 6) ? new(elastic_matrix) : error("Size should be 6x6!")
end

compliance_matrix(c::CrystalSystem) = inv(c.elastic_matrix)

stability_conditions(c::CrystalSystem) = "No stability conditions defined for `$(typeof(c))!`"
stability_conditions(c::CubicSystem) = [
    "C_{11} > | C_{12} |",
    "C_{11} + 2 C_{12} > 0",
    "C_{44} > 0"
]
stability_conditions(c::HexagonalSystem) = [
    "C_{11} > | C_{12} |",
    "2 C_{13}^2 < C_{33} (C_{11} + C_{12})",
    "C_{44} > 0",
    "C_{66} > 0"
]


end