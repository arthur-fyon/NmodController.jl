# Simulating conductance based models
## Simulating uncontrolled conductance based models
The most intersting thing about conductance based models kicks in when integrating their ODEs. Indeed, the voltage solution can reproduce neuronal activity and is widely used in either computational neuroscience or neuromorphic engineering. *NmodController.jl* can be used to simulate conductance based models coupled with [*DifferentialEquations.jl*](https://github.com/SciML/DifferentialEquations.jl). This is a suite for numerically solving differential equations written in Julia. The purpose of this package is to supply efficient Julia implementations of solvers for various differential equations.

The way *NmodController.jl* helps with simulating conductance based model is that it provides an easy way to write the ODE function necessary for *DifferentialEquations.jl* solvers with the function `writeUncontrolledODEs()`. As arguments, you only have to specify your model. Note that the ODE function will be written in a newly created *.jl* file that you can specify with the optional argument `filename`.

### Example 1
The next few lines of code show how to write the ODE function file using a predefined conductance based model `neuron`.

```julia
using NmodController
writeUncontrolledODEs(neuron, filename="neuron_ODE.jl")
```

To define your ODE problem, you have to include your newly created *.jl* file and to put the ODE function, which is by default `CB_ODE`, in argument of the `ODEproblem()` function from the *DifferentialEquations.jl* suite along with the simulation parameters, initial conditions and time span over which you want to integrate. Finally, just solve with `solve()` from the *DifferentialEquations.jl* suite.

## Simulating controlled conductance based models
In Fyon et al., 2023 "Reliable neuromodulation from adaptive control of ion channel expression" is proposed a first robust neuromodulation adaptive controller that can be used with any conductance based model. In brief, it allows you to control $p$ DICs at any voltage (here, we only use threshold voltage since DICs at such voltage shape neuronal activity) with $n$ ionic currents. Thus, by adequately choosing the neuromodulated maximum ion channel conductances and the DICs you want to control, any physiological firing pattern can be achieved by inputting appropriate values of the DICs. Note that $p \leq n$ in order to be the linear system behind the controller to have a solution.

The way *NmodController.jl* helps with simulating neuromodulated conductance based model is that it provides an easy way to write the ODE function necessary for *DifferentialEquations.jl* solvers with the function `writeControlledODEs()`. As arguments, you only have to specify your model, the names of the conductance you want to neuromodulate and the names of the DICs you want to control. Note that the ODE function will be written in a newly created *.jl* file that you can specify with the optional argument `filename`.

### Example 2
The next few lines of code show how to write the ODE function file using a predefined conductance based model `neuron`.

```julia
using NmodController
writeControlledODEs(neuron, ["ion1", "ion2", "ion3"], ["f", "s", "u"], filename="ControlledNeuron_ODE.jl")
```

To define your ODE problem, you have to include your newly created *.jl* file and to put the ODE function, which is by default `Controlled_CB_ODE`, in argument of the `ODEproblem()` function from the *DifferentialEquations.jl* suite along with the simulation parameters, initial conditions and time span over which you want to integrate. Finally, just solve with `solve()` from the *DifferentialEquations.jl* suite.