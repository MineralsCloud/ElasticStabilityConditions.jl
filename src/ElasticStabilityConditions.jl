module ElasticStabilityConditions

using LinearAlgebra: diag

using Crystallography
using StaticArrays: FieldVector, FieldMatrix, FieldArray, SVector, SHermitianCompact

export TensorStress,
    TensorStrain,
    TensorStiffness,
    TensorCompliance,
    EngineeringStress,
    EngineeringStrain,
    EngineeringCompliance,
    EngineeringStiffness
export stability_conditions, isstable

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

struct TensorStiffness{T} <: FieldArray{Tuple{3,3,3,3},T,4}
    xxxx::T
    yxxx::T
    zxxx::T
    xyxx::T
    yyxx::T
    zyxx::T
    xzxx::T
    yzxx::T
    zzxx::T
    xxyx::T
    yxyx::T
    zxyx::T
    xyyx::T
    yyyx::T
    zyyx::T
    xzyx::T
    yzyx::T
    zzyx::T
    xxzx::T
    yxzx::T
    zxzx::T
    xyzx::T
    yyzx::T
    zyzx::T
    xzzx::T
    yzzx::T
    zzzx::T
    xxxy::T
    yxxy::T
    zxxy::T
    xyxy::T
    yyxy::T
    zyxy::T
    xzxy::T
    yzxy::T
    zzxy::T
    xxyy::T
    yxyy::T
    zxyy::T
    xyyy::T
    yyyy::T
    zyyy::T
    xzyy::T
    yzyy::T
    zzyy::T
    xxzy::T
    yxzy::T
    zxzy::T
    xyzy::T
    yyzy::T
    zyzy::T
    xzzy::T
    yzzy::T
    zzzy::T
    xxxz::T
    yxxz::T
    zxxz::T
    xyxz::T
    yyxz::T
    zyxz::T
    xzxz::T
    yzxz::T
    zzxz::T
    xxyz::T
    yxyz::T
    zxyz::T
    xyyz::T
    yyyz::T
    zyyz::T
    xzyz::T
    yzyz::T
    zzyz::T
    xxzz::T
    yxzz::T
    zxzz::T
    xyzz::T
    yyzz::T
    zyzz::T
    xzzz::T
    yzzz::T
    zzzz::T
end

struct TensorCompliance{T} <: FieldArray{Tuple{3,3,3,3},T,4}
    xxxx::T
    yxxx::T
    zxxx::T
    xyxx::T
    yyxx::T
    zyxx::T
    xzxx::T
    yzxx::T
    zzxx::T
    xxyx::T
    yxyx::T
    zxyx::T
    xyyx::T
    yyyx::T
    zyyx::T
    xzyx::T
    yzyx::T
    zzyx::T
    xxzx::T
    yxzx::T
    zxzx::T
    xyzx::T
    yyzx::T
    zyzx::T
    xzzx::T
    yzzx::T
    zzzx::T
    xxxy::T
    yxxy::T
    zxxy::T
    xyxy::T
    yyxy::T
    zyxy::T
    xzxy::T
    yzxy::T
    zzxy::T
    xxyy::T
    yxyy::T
    zxyy::T
    xyyy::T
    yyyy::T
    zyyy::T
    xzyy::T
    yzyy::T
    zzyy::T
    xxzy::T
    yxzy::T
    zxzy::T
    xyzy::T
    yyzy::T
    zyzy::T
    xzzy::T
    yzzy::T
    zzzy::T
    xxxz::T
    yxxz::T
    zxxz::T
    xyxz::T
    yyxz::T
    zyxz::T
    xzxz::T
    yzxz::T
    zzxz::T
    xxyz::T
    yxyz::T
    zxyz::T
    xyyz::T
    yyyz::T
    zyyz::T
    xzyz::T
    yzyz::T
    zzyz::T
    xxzz::T
    yxzz::T
    zxzz::T
    xyzz::T
    yyzz::T
    zyzz::T
    xzzz::T
    yzzz::T
    zzzz::T
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

struct EngineeringCompliance{T} <: FieldMatrix{6,6,T}
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

struct EngineeringStiffness{T} <: FieldMatrix{6,6,T}
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

stability_conditions(::Cubic) =
    ["C_{11} > | C_{12} |", "C_{11} + 2 C_{12} > 0", "C_{44} > 0"]
stability_conditions(::Hexagonal) = [
    "C_{11} > | C_{12} |",
    "2 C_{13}^2 < C_{33} (C_{11} + C_{12})",
    "C_{44} > 0",
    "C_{66} > 0",
]

function isstable(::Cubic, c::EngineeringStiffness)
    c11, c12, c44 = c[1, 1], c[1, 2], c[1, 4]
    return all([  # Must satisfy all criteria!
        c11 > abs(c12),
        c11 + 2 * c12 > 0,
        c44 > 0,
    ])
end # function isstable
function isstable(::Hexagonal, c::EngineeringStiffness)
    c11, c12, c13, c33, c44, c66 = c[1, 1], c[1, 2], c[1, 3], c[3, 3], c[4, 4], c[6, 6]
    return all([  # Must satisfy all criteria!
        c11 > abs(c12),
        2 * c13^2 < c33 * (c11 + c12),
        c44 > 0,
        c66 > 0,
    ])
end # function isstable
function isstable(::Tetragonal, c::EngineeringStiffness)
    c11, c12, c13, c16, c33, c44, c66 =
        c[0, 0], c[0, 1], c[0, 2], c[0, 5], c[2, 2], c[3, 3], c[5, 5]
    if c16 == 0  # Tetragonal (I) class
        return [c11 > abs(c12), 2 * c13^2 < c33 * (c11 + c12), c44 > 0, c66 > 0]
    end
    return [  # Tetragonal (II) class
        c11 > abs(c12),
        2 * c13^2 < c33 * (c11 + c12),
        c44 > 0,
        2 * c16^2 < c66 * (c11 - c12),
    ]
end # function isstable
function isstable(::Trigonal, c::EngineeringStiffness)
    c11, c12, c13, c14, c15, c33, c44, c66 =
        c[0, 0], c[0, 1], c[0, 2], c[0, 3], c[0, 4], c[2, 2], c[3, 3], c[5, 5]
    if c15 == 0  # Rhombohedral (I) class
        return [
            c11 > abs(c12),
            c44 > 0,
            c13^2 < 0.5 * c33 * (c11 + c12),
            c14^2 < 0.5 * c44 * (c11 - c12),
            0.5 * c44 * (c11 - c12) == c44 * c66,
        ]
    end
    return [  # Rhombohedral (II) class
        c11 > abs(c12),
        c44 > 0,
        c13^2 < 0.5 * c33 * (c11 + c12),
        c14^2 + c15^2 < 0.5 * c44 * (c11 - c12),
        0.5 * c44 * (c11 - c12) == c44 * c66,
    ]
end # function isstable
function isstable(::Orthorhombic, c::EngineeringStiffness)
    c11, c22, c33, c44, c55, c66 = diag(c)
    c12, c13, c23 = c[0, 1], c[0, 2], c[1, 2]
    return [
        c11 > 0,
        c11 * c22 > c12^2,
        c11 * c22 * c33 + 2 * c12 * c13 * c23 > c11 * c23^2 + c22 * c13^2 + c33 * c12^2,
        c44 > 0,
        c55 > 0,
        c66 > 0,
    ]
end # function isstable
function isstable(::Monoclinic, c::EngineeringStiffness)
    c11, c22, c33, c44, c55, c66 = diag(c)
    c12, c13, c15, c23, c25, c35, c46 =
        c[0, 1], c[0, 2], c[0, 4], c[1, 2], c[1, 4], c[2, 4], c[3, 5]
    g =
        c11 * c22 * c33 - c11 * c23 * c23 - c22 * c13 * c13 - c33 * c12 * c12 +
        2 * c12 * c13 * c23
    return [
        all(_ > 0 for _ in (c11, c22, c33, c44, c55, c66)),
        c11 + c22 + c33 + 2 * (c12 + c13 + c23) > 0,
        c33 * c55 - c35^2 > 0,
        c44 * c66 - c46^2 > 0,
        c22 + c33 - 2 * c23 > 0,
        c22 * (c33 * c55 - c35^2) + 2 * c23 * c25 * c35 - c23^2 * c55 - c25^2 * c33 > 0,
        2 * (
            c15 * c25 * (c33 * c12 - c13 * c23) +
            c15 * c35 * (c22 * c13 - c12 * c23) +
            c25 * c35 * (c11 * c23 - c12 * c13)
        ) - (
            c15 * c15 * (c22 * c33 - c23^2) +
            c25 * c25 * (c11 * c33 - c13^2) +
            c35 * c35 * (c11 * c22 - c12^2)
        ) + c55 * g > 0,
    ]
end # function isstable
isstable(C::CrystalSystem, s::EngineeringCompliance) = isstable(C, inv(s))

function Base.convert(::Type{<:TensorStress}, s::EngineeringStress)
    return TensorStress(SHermitianCompact(SVector(s.xx, s.xy, s.xz, s.yy, s.yz, s.zz)))
end # function Base.convert
function Base.convert(::Type{<:EngineeringStress}, s::TensorStress)
    return EngineeringStress(s.xx, s.yy, s.zz, s.yz, s.xz, s.xy)
end # function Base.convert
function Base.convert(::Type{<:TensorStrain}, e::EngineeringStrain)
    return TensorStrain(SHermitianCompact(SVector(
        e.xx,
        e.xy / 2,
        e.xz / 2,
        e.yy,
        e.yz / 2,
        e.zz,
    )))
end # function Base.convert
function Base.convert(::Type{<:EngineeringStrain}, e::TensorStrain)
    return EngineeringStrain(e.xx, e.yy, e.zz, 2 * e.yz, 2 * e.xz, 2 * e.xy)
end # function Base.convert

Base.inv(c::EngineeringStiffness) = EngineeringCompliance(inv(collect(c)))
Base.inv(s::EngineeringCompliance) = EngineeringStiffness(inv(collect(s)))

end
