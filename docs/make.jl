using PioneerBenchmarks
using Documenter

DocMeta.setdocmeta!(PioneerBenchmarks, :DocTestSetup, :(using PioneerBenchmarks); recursive=true)

makedocs(;
    modules=[PioneerBenchmarks],
    authors="Nathan Wamsley",
    repo="https://github.com/nwamsley1/PioneerBenchmarks.jl/blob/{commit}{path}#{line}",
    sitename="PioneerBenchmarks.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://nwamsley1.github.io/PioneerBenchmarks.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/nwamsley1/PioneerBenchmarks.jl",
)
