using Documenter

makedocs(
    sitename = "NmodController",
    format = Documenter.HTML(),
    modules = [NmodController]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/arthur-fyon/NmodController.jl"
)
