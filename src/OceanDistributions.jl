module OceanDistributions

using CSV, DataFrames
export readoceandistribution

"""
    readoceandistribution(fil::String)

Read a distribution from csv file

```
Tcensus=readoceandistribution("examples/Tcensus.txt")
```
"""
function readoceandistribution(file::String)
    return CSV.read(file,DataFrame)
end

end # module
