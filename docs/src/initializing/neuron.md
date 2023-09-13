# Initializing a conductance based model

Once all the ionic current data structures have been initialized, any complete conductance based model can be descrbied. The voltage equation for such model, without any external applied current, writes

$\dot{V} = (1/C) \cdot (-\sum_\mathrm{ion} I_\mathrm{ion} - I_\mathrm{leak})$

where $C$ is the membrane capacitance and $I_\mathrm{leak}$ writes

$I_\mathrm{leak} = g_\mathrm{leak} (V - E_\mathrm{leak}).$

All the other equations in the model consists in the gating variables dynamics previously described.

In *NmodController.jl*, a conductance based model can be contained in a `NeuronCB` type. To help initializing such data structure, calling `initializeNeuronModel()` with appropriate arguments is strongly recommended.

### Example 1

The next few lines of code show how to initialize a conductance based model with two ionic currents: a fast sodium and a rectified delayed potassium, i.e. the original Hodgkin and Huxley model. This voltage equation writes

$\dot{V} = (1/C) \cdot (-I_\mathrm{Na} -  I_\mathrm{Kd} - I_\mathrm{leak}).$

```julia
# First wrapping all ionic currents in a vector
ionCurrents = [NaCurrent, KdCurrent]

# Initializing leakage reversal potential and leakage conductance
Eleak = -50.
gleak = 0.01

# Initializing all the maximum ion channel conductances
bar_g = [100., 10.]

using NmodController
HHmodel = initializeNeuronModel(ionCurrents, C=0.1, leakageConductance=gleak, reversaleLeakagePotential=Eleak, maximumConductances=bar_g)
```

Note that the argument `maximumConductances` is optional and correspond to all the maximum ion channel conductances $\bar{g}_\mathrm{ion}$ wrapped in a vector. `maximumConductances` will be filled with `NaN` if not provided in input.

## Intracellular calcium dynamic in a conductance based model