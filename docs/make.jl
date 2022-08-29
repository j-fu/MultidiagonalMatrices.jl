using Documenter, MultidiagonalMatrices, StaticArrays, LinearAlgebra, SparseArrays

function mkdocs()
    makedocs(sitename="MultidiagonalMatrices",
             modules = [MultidiagonalMatrices],
             doctest = true,
             clean = false,
             authors = "J. Fuhrmann",
             repo="https://github.com/j-fu/MultidiagonalMatrices.jl",
             pages=[
                 "Home"=>"index.md",
                 "API" => "api.md",
                 "Internal" => "internal.md"
             ])
end

mkdocs()

deploydocs(repo = "github.com/j-fu/MultidiagonalMatrices.jl.git",branch="main")
