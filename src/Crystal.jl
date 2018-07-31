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

abstract type CrystalSystem end
abstract type CubicSystem <: CrystalSystem end
abstract type HexagonalSystem <: CrystalSystem end

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

satisfy_stability_conditions(c::CrystalSystem) = false
function satisfy_stability_conditions(cub::CubicSystem)
    c = cub.elastic_matrix
    c11, c12, c44 = c[1, 1], c[1, 2], c[1, 4]
    criteria = [
        c11 > abs(c12),
        c11 + 2 * c12 > 0,
        c44 > 0
    ]
    all(criteria)
end
function satisfy_stability_conditions(hex::HexagonalSystem)
    c = hex.elastic_matrix
    c11, c12, c13, c33, c44, c66 = c[1, 1], c[1, 2], c[1, 3], c[3, 3], c[4, 4], c[6, 6]
    criteria = [
        c11 > abs(c12),
        2 * c13^2 < c33 * (c11 + c12),
        c44 > 0,
        c66 > 0
    ]
    all(criteria)
end

end
