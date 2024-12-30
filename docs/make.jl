using Documenter

makedocs(
    sitename = "VectorPSFs.jl",
    pages = [
        "Home" => "index.md",
        "References" => "apiref.md"
    ],
    modules = [VectorPSFs],
    format = Documenter.HTML()
)
deploydocs(
    repo = "github.com/IvanKuznetsoff/VectorPSFs.jl.git",
)