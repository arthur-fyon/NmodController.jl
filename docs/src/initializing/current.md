# Initializing an ionic current

In conductance based models, an ionic current $I_\mathrm{ion}$ is described by its voltage dependent conductance $g_\mathrm{ion}$ and Nernst reversal potential $E_\mathrm{ion}$

$I_\mathrm{ion} = g_\mathrm{ion}(V) \cdot (V - E_\mathrm{ion}).$

This voltage dependent conductance is defined by its maximum value $\bar{g}_\mathrm{ion}$ and at most two gating variables, one activation $m_\mathrm{ion}(V)$ and one inactivation $h_\mathrm{ion}(V)$, that varies between 0 and 1. In the case of a current with two gating variables, this writes

$g_\mathrm{ion}(V) = \bar{g}_\mathrm{ion} \cdot m^{a}_\mathrm{ion}(V) \cdot h^{b}_\mathrm{ion}(V).$

Note that the maximum ion channel conductance of the current will not be stocked inside the current data structure, but rather in the conductance based model data structure.

Each gating variable dynamic follow a basic first order ODE where both the time constant $\tau_{m_\mathrm{ion}}$ and the converging value $m_{\mathrm{ion}_\infty}$ are voltage dependent: 

$\dot{m}_\mathrm{ion} = \frac{m_{\mathrm{ion}_\infty}(V) - m_\mathrm{ion}}{\tau_{m_\mathrm{ion}}(V)}.$

In *NmodController.jl*, an ionic current can be contained in a `IonCurrent` type. To help initializing such data structure, calling `initializeCurrent()` with appropriate arguments is strongly recommended.

### Example 1

The next few lines of code show how to initialize a sodium current with two gating variables. This current writes

$I_\mathrm{Na} = \bar{g}_\mathrm{Na} \cdot m^{3}_\mathrm{Na} \cdot h_\mathrm{Na} \cdot (V - E_\mathrm{Na}).$

```@example
# First initializing the converging values and time constants functions
boltz(V, A, B) = 1 / (1 + exp((V+A) / B))
tauX(V, A, B, D, E) = A - B / (1 + exp((V+D) / E))

mNa_inf(V) = boltz(V, 25.5, -5.29)
tau_mNa(V) = tauX(V, 1.32, 1.26, 120., -25.)
hNa_inf(V) = boltz(V, 48.9, 5.18)
tau_hNa(V) = (0.67 / (1 + exp((V+62.9)/-10.0))) * (1.5 + 1 / (1+exp((V+34.9)/3.6)))

# Initializing Nernst reversal potential
ENa = 50.

using NmodController
NaCurrent = initializeCurrent("Na", ENa, numberOfGatings=2, exponents=[3, 1],
    activationSteadyStateGating=mNa_inf, activationTimeConstant=tau_mNa,
    inactivationSteadyStateGating=hNa_inf, inactivationTimeConstant=tau_hNa)
```

### Example 2

Sometimes, the gating variable converging values may depend on the intracellular calcium. The next few lines of code show how to initialize a calcium controlled potassium current with one gating variable. This current writes

$I_\mathrm{KCa} = \bar{g}_\mathrm{KCa} \cdot m^{4}_\mathrm{KCa} \cdot (V - E_\mathrm{KCa}).$

Where the dynamic of $m_\mathrm{KCa}$ is described by

$\dot{m}_\mathrm{KCa} = \frac{m_{\mathrm{KCa}_\infty}(V, Ca) - m_\mathrm{KCa}}{\tau_{m_\mathrm{KCa}}(V)}.$

As you can see, $m_{\mathrm{KCa}_\infty}$ depends on both the voltage and the calcium, this must be notified to `initializeCurrent()`.

```julia
# First initializing the converging value and time constant functions
boltz(V, A, B) = 1 / (1 + exp((V+A) / B))
tauX(V, A, B, D, E) = A - B / (1 + exp((V+D) / E))

mKCa_inf(V, Ca) = (Ca / (Ca+3)) * (1 / (1+exp((V+28.3)/-12.6)))
tau_mKCa(V) = tauX(V, 90.3, 75.1, 46., -22.7)

# Initializing Nernst reversal potential
EK = -80.

using NmodController
KCaCurrent = initializeCurrent("KCa", EK, exponents=4,
    activationSteadyStateGating=mKCa_inf, activationTimeConstant=tau_mKCa,
    calciumDependency=true)
```

### Example 3

Another specific case is where the time constant may not depend on the voltage. In such case, just provide the time constant as a `Float64` or a `Int64` in argument of `initializeCurrent()`. Moreover, the converging value of the gating variable may be magnesium dependent. In such case, do as the following. The next few lines of code show how to initialize an instantaneous magnesium dependent NMDA current with one gating variable. This current writes

$I_\mathrm{NMDA} = \bar{g}_\mathrm{NMDA} \cdot m_\mathrm{NMDA} \cdot (V - E_\mathrm{NMDA}).$

Where the dynamic of $m_\mathrm{NMDA}$ is described by

$m_\mathrm{NMDA} = m_{\mathrm{NMDA}_\infty}(V, Mg).$

As you can see, $m_{\mathrm{NMDA}_\infty}$ depends on both the voltage and the magnesium, this must be notified to `initializeCurrent()`

```julia
# First initializing the converging value and time constant functions
NMDA_inf(V, Mg) = 1 / (1 + Mg*exp(-0.08*V)/10.)
tau_NMDA(V) = 1e-10 # Never equal to zero, this might be fixed later on

# Initializing Nernst reversal potential
ENMDA = 0.

using NmodController
NMDACurrent = initializeCurrent("NMDA", ENMDA, exponents=1,
    activationSteadyStateGating=NMDA_inf, activationTimeConstant=tau_NMDA,
    MgDependency=true)
```