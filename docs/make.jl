Pkg.activate("./docs");

using Documenter, DistributionsLH, FilesLH

makedocs(
    modules = [DistributionsLH],
    format = Documenter.HTML(; prettyurls = get(ENV, "CI", nothing) == "true"),
    authors = "hendri54",
    sitename = "DistributionsLH.jl",
    pages = Any["index.md"]
    # strict = true,
    # clean = true,
    # checkdocs = :exports,
)

pkgDir = rstrip(normpath(@__DIR__, ".."), '/');
@assert endswith(pkgDir, "DistributionsLH")
deploy_docs(pkgDir);

Pkg.activate(".");


# deploydocs(
#     repo = "github.com/hendri54/DistributionsLH.jl.git",
#     push_preview = true
# )
