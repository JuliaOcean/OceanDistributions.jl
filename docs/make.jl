using Documenter, OceanDistributions

makedocs(;
    modules=[OceanDistributions],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/gaelforget/OceanDistributions.jl/blob/{commit}{path}#L{line}",
    sitename="OceanDistributions.jl",
    authors="gaelforget <gforget@mit.edu>",
    assets=String[],
)

deploydocs(;
    repo="github.com/gaelforget/OceanDistributions.jl",
)
