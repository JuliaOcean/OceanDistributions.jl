using Documenter, OceanDistributions
import PlutoSliderServer

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

lst=("one_dim_diffusion.jl",)
for i in lst
    fil_in=joinpath(@__DIR__,"..", "examples",i)
    fil_out=joinpath(@__DIR__,"build", i[1:end-2]*"html")
    PlutoSliderServer.export_notebook(fil_in)
    mv(fil_in[1:end-2]*"html",fil_out)
    cp(fil_in,fil_out[1:end-4]*"jl")
end

deploydocs(;
    repo="github.com/gaelforget/OceanDistributions.jl",
)
