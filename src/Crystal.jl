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
end

function compliance_matrix(m::CrystalSystem)
    return inv(m)
end



end