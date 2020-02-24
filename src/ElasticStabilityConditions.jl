module ElasticStabilityConditions

using Crystallography
using StaticArrays: FieldVector, FieldMatrix, FieldArray, SVector, SHermitianCompact

export TensorStress,
    TensorStrain,
    TensorStiffness,
    TensorCompliance,
    EngineeringStress,
    EngineeringStrain,
    Compliance,
    Stiffness
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

stability_conditions(::Cubic) =
    ["C_{11} > | C_{12} |", "C_{11} + 2 C_{12} > 0", "C_{44} > 0"]
stability_conditions(::Hexagonal) = [
    "C_{11} > | C_{12} |",
    "2 C_{13}^2 < C_{33} (C_{11} + C_{12})",
    "C_{44} > 0",
    "C_{66} > 0",
]

function isstable(::Cubic, c::Stiffness)
    c11, c12, c44 = c[1, 1], c[1, 2], c[1, 4]
    return all([  # Must satisfy all criteria!
        c11 > abs(c12),
        c11 + 2 * c12 > 0,
        c44 > 0,
    ])
end
function isstable(::Hexagonal, c::Stiffness)
    c11, c12, c13, c33, c44, c66 = c[1, 1], c[1, 2], c[1, 3], c[3, 3], c[4, 4], c[6, 6]
    return all([  # Must satisfy all criteria!
        c11 > abs(c12),
        2 * c13^2 < c33 * (c11 + c12),
        c44 > 0,
        c66 > 0,
    ])
end
isstable(C::CrystalSystem, s::Compliance) = isstable(C, inv(s))

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

Base.inv(c::Stiffness) = Compliance(inv(collect(c)))
Base.inv(s::Compliance) = Stiffness(inv(collect(s)))

end
