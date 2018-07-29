"""
# module Crystal.jl

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
    CrystalSystem(elastic_matrix) = size(elastic_matrix) == (2, 2) ? new(elastic_matrix) : error("Size should be 6x6!")
end

function compliance_matrix(m::CrystalSystem)
    return inv(m.elastic_matrix)
end

println(compliance_matrix(CrystalSystem([1 2;3 4])))

end