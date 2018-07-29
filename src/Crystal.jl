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

compliance_matrix(m::CrystalSystem) = inv(m.elastic_matrix)

end