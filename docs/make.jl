using PDEfind
using Documenter

DocMeta.setdocmeta!(PDEfind, :DocTestSetup, :(using PDEfind); recursive=true)

makedocs(;
    modules=[PDEfind],
    authors="Natan Dominko Kobilica <natandominko@gmail.com>",
    sitename="PDEfind.jl",
    format=Documenter.HTML(;
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
