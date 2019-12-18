using ElasticStabilityConditions
using Documenter

makedocs(;
    modules=[ElasticStabilityConditions],
    authors="Qi Zhang <singularitti@outlook.com>",
    repo="https://github.com/MineralsCloud/ElasticStabilityConditions.jl/blob/{commit}{path}#L{line}",
    sitename="ElasticStabilityConditions.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://MineralsCloud.github.io/ElasticStabilityConditions.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/MineralsCloud/ElasticStabilityConditions.jl",
)
