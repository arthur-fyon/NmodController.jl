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

When at least one ionic current is calcium dependent, an additional ODE has to be added to the conductance based model to describe intracellular calcium dynamics. Such equation generally writes

$\tau_{Ca} \cdot \dot{\[Ca\]} = \sum_{\mathrm{ion}\,Ca} e_{\mathrm{ion}\,Ca} I_{\mathrm{ion}\,Ca} - \[Ca\] + Ca_\mathrm{eq}$

where $I_{\mathrm{ion}\,Ca}$ is a calcium current of the model and $e_{\mathrm{ion}\,Ca}$ is its associated coefficient.

In *NmodController.jl*, such intracellular calcium dynamics can be contained in a `CalciumDynamic` type. To help initializing such data structure, calling `initializeCalciumDynamics()` with appropriate arguments is strongly recommended.

### Example 2
The next few lines of code show how to initialize a conductance based model with five ionic currents: a fast sodium, a rectified delayed potassium, a calcium controlled potassium current, a T-type calcium current and a slow calcium current. Note that the calcium controlled potassium current is a special current that depends on both the voltage and intracellular calcium. This voltage equation writes

$\dot{V} = (1/C) \cdot (-I_\mathrm{Na} -  I_\mathrm{Kd} -  I_\mathrm{KCa} -  I_\mathrm{CaS} -  I_\mathrm{CaT} - I_\mathrm{leak}).$

The associated calcium dynamic writes

$\tau_{Ca} \cdot \dot{\[Ca\]} = -0.94\cdot  I_{CaT} -0.94\cdot  I_{CaS} - \[Ca\] + Ca_\mathrm{eq}$

with $\tau_{Ca} = 20$ and $Ca_\mathrm{eq} = 0.05$.

```julia
# First wrapping all ionic currents in a vector
ionCurrents = [NaCurrent, KdCurrent, KCaCurrent, CaTCurrent, CaSCurrent]

# Initializing leakage reversal potential and leakage conductance
Eleak = -50.
gleak = 0.01

# Initializing all the maximum ion channel conductances
bar_g = [100., 10., 10., 1., 1.]

# Initialize the current dynamic data structure
CaDyn = initializeCalciumDynamics(["CaT", "CaS"], [-0.94, -0.94], 0.05, 20)

using NmodController
HHmodel = initializeNeuronModel(ionCurrents, C=0.1, calciumDynamics=CaDyn, leakageConductance=gleak, reversaleLeakagePotential=Eleak, maximumConductances=bar_g)
```