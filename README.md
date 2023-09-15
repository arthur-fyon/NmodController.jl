[![Docs](https://img.shields.io/badge/docs-NmodController-blue.svg)](https://arthur-fyon.github.io/NmodController.jl/)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/arthur-fyon/NmodController.jl?svg=true)](https://ci.appveyor.com/project/arthur-fyon/NmodController-jl)
[![Coverage](https://codecov.io/gh/arthur-fyon/NmodController.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/arthur-fyon/NmodController.jl)

# NmodController.jl: A simple way to use conductance based models and to simulate them

This is a suite for numerically simulating conductance based model when coupled with [*DifferentialEquations.jl*](https://github.com/SciML/DifferentialEquations.jl). The purpose of this package is to supply easy definition and simulation of any conductance based model. Moreover, it permits easy broadcasting of Dynamics Input Conductance theory ([paper](https://www.eneuro.org/content/2/1/ENEURO.0031-14.2015)) and of a robust neuromodulation controller ([paper](https://orbi.uliege.be/handle/2268/304427)).

## Installation

### Users
1) Download [Julia v1.6.X](https://julialang.org/downloads/) or later, if you haven't already.
1) Add the NmodController module entering the following at the Julia REPL `]add https://github.com/arthur-fyon/NmodController.jl`.

### Developers
1) Clone the NmodController module to `username/.julia/dev/`.
2) Enter the package manager in REPL by pressing `]` then add the package by typing `dev NmodController` rather than `add NmodController`.

## [Documentation](https://arthur-fyon.github.io/NmodController.jl/)