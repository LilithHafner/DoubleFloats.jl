using Documenter, DoubleFloats

makedocs(
    modules = [DoubleFloats],
    sitename = "DoubleFloats.jl",
    authors = "Jeffrey Sarnoff, Simon Byrne, Sascha Timme, and other contributors",
    pages = Any[
        "Overview" => "index.md",
        "Construction" => "construction.md",
        "Characteristics" => "characteristics.md",
        "String, Show, Parse" => "stringshowparse.md",
        "Capabilities" => "capabilities.md",
        "Additional Capabilities" => "special.md",
        "Linear Algebra" => "linearalgebra.md",
        "Appropriate Uses" => "appropriate.md",
        "References" => "references.md"
    ]
)

deploydocs(
    repo = "github.com/JuliaMath/DoubleFloats.jl.git",
    target = "build"
)
