using Documenter, OceanWaterMasses

makedocs(;
    modules=[OceanWaterMasses],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/gaelforget/OceanWaterMasses.jl/blob/{commit}{path}#L{line}",
    sitename="OceanWaterMasses.jl",
    authors="gaelforget <gforget@mit.edu>",
    assets=String[],
)

deploydocs(;
    repo="github.com/gaelforget/OceanWaterMasses.jl",
)
