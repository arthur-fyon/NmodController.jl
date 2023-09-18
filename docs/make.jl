using Documenter

using SymPy

using NmodController

makedocs(
    sitename = "NmodController",
    format = Documenter.HTML(),
    pages=[
        "Home" => "index.md",
        "Initializing a conductance based model" =>
                    ["How to initialize an ionic current?"          => "initializing/current.md",
                     "How to initialize a conductance based model?" => "initializing/neuron.md",
                     "Examples of existing models"                  => "initializing/examples.md",],
        "Computing dynamical input conductances (DICs)" =>
                    ["How to compute DICs and sensitivity matrix?" => "computing/DICs.md",
                     "Examples of existing models"                 => "computing/examples.md",],
        "Simulating with DifferentialEquations.jl" =>
                    ["How to write ODE functions file?" => "simulating/files.md",
                     "Examples of existing models"      => "simulating/examples.md",],
        "Types/Methods/Functions" => "typesMethodsFunctions.md",
    ],
    modules = Module[NmodController]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/arthur-fyon/NmodController.jl"
)
