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

using Crystallography
using StaticArrays: FieldVector, FieldMatrix, FieldArray

struct TensorStress{T} <: FieldMatrix{3,3,T}
    xx::T
    yx::T
    zx::T
    xy::T
    yy::T
    zy::T
    xz::T
    yz::T
    zz::T
end

struct TensorStrain{T} <: FieldMatrix{3,3,T}
    xx::T
    yx::T
    zx::T
    xy::T
    yy::T
    zy::T
    xz::T
    yz::T
    zz::T
end

struct EngineeringStress{T} <: FieldVector{6,T}
    xx::T
    yy::T
    zz::T
    yz::T
    xz::T
    xy::T
end

struct EngineeringStrain{T} <: FieldVector{6,T}
    xx::T
    yy::T
    zz::T
    yz::T
    xz::T
    xy::T
end

struct Compliance{T} <: FieldMatrix{6,6,T}
    xxxx::T
    yyxx::T
    zzxx::T
    yzxx::T
    xzxx::T
    xyxx::T
    xxyy::T
    yyyy::T
    zzyy::T
    yzyy::T
    xzyy::T
    xyyy::T
    xxzz::T
    yyzz::T
    zzzz::T
    yzzz::T
    xzzz::T
    xyzz::T
    xxyz::T
    yyyz::T
    zzyz::T
    yzyz::T
    xzyz::T
    xyyz::T
    xxxz::T
    yyxz::T
    zzxz::T
    yzxz::T
    xzxz::T
    xyxz::T
    xxxy::T
    yyxy::T
    zzxy::T
    yzxy::T
    xzxy::T
    xyxy::T
end

struct Stiffness{T} <: FieldMatrix{6,6,T}
    xxxx::T
    yyxx::T
    zzxx::T
    yzxx::T
    xzxx::T
    xyxx::T
    xxyy::T
    yyyy::T
    zzyy::T
    yzyy::T
    xzyy::T
    xyyy::T
    xxzz::T
    yyzz::T
    zzzz::T
    yzzz::T
    xzzz::T
    xyzz::T
    xxyz::T
    yyyz::T
    zzyz::T
    yzyz::T
    xzyz::T
    xyyz::T
    xxxz::T
    yyxz::T
    zzxz::T
    yzxz::T
    xzxz::T
    xyxz::T
    xxxy::T
    yyxy::T
    zzxy::T
    yzxy::T
    xzxy::T
    xyxy::T
end

stability_conditions(c::CrystalSystem) = "No stability conditions defined for `$(typeof(c))!`"
stability_conditions(::Cubic) = [
    "C_{11} > | C_{12} |",
    "C_{11} + 2 C_{12} > 0",
    "C_{44} > 0"
]
stability_conditions(::Hexagonal) = [
    "C_{11} > | C_{12} |",
    "2 C_{13}^2 < C_{33} (C_{11} + C_{12})",
    "C_{44} > 0",
    "C_{66} > 0"
]

satisfy_stability_conditions(::CrystalSystem) = false
function satisfy_stability_conditions(cub::Cubic)
    c = cub.elastic_matrix
    c11, c12, c44 = c[1, 1], c[1, 2], c[1, 4]
    criteria = [
        c11 > abs(c12),
        c11 + 2 * c12 > 0,
        c44 > 0
    ]
    all(criteria)
end
function satisfy_stability_conditions(hex::Hexagonal)
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
