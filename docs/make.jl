using Documenter
using DS9ui

makedocs(
    sitename="DS9ui.jl",
    authors="Marco Lombardi",
    format=Documenter.HTML(
        #assets=["assets/favicon/favicon.ico"]
    ),
    pages=[
        "Home" => "index.md",
        "API" => Any[
            "Connection"=>"connection.md",
            "Interface"=>"interface.md",
            "Analysis"=>"analysis.md",
        ],
    ])

deploydocs(
    repo = "github.com/astrozot/DS9ui.jl",
    versions = nothing
)