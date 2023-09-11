using Documenter
using NmodController

makedocs(
    sitename = "NmodController",
    format = Documenter.HTML(),
    modules = [NmodController]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
