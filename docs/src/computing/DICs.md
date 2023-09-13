# Computations on conductance based model
## Computing dynamic input conductances (DICs)
DICs are three voltage-dependent conductance curves $g_\mathrm{f}(V),g_\mathrm{s}(V),g_\mathrm{u}(V)$ that can  be computed as linear functions of the maximal conductance vector $\bar{g}_\mathrm{ion}$ of the neuron model at each $V$. By wrapping the three dynamic input conductances in a vector $g_\mathrm{DIC}(V)$, the linear combination writes

$g_\mathrm{DIC}(V) = S(V) \cdot \bar{g}_\mathrm{ion}$

where $S(V)$ is a sensitivity matrix. These mathematical function can be built using Drion et al., 2015 "Dynamic input conductances shape neuronal spiking".

Because of the specific feedback structure of conductance-based models, DICs shape neuronal spiking behavior and the three DICs differ in the timescale at which this shaping happens: fast $g_\mathrm{f}(V)$, slow $g_\mathrm{s}(V)$, and ultraslow $g_\mathrm{u}(V)$. Values and signs of the DICs at specific voltages, mainly the threshold voltage $V_\mathrm{th}$, reliably determine the neuronal firing pattern. For instance, a negative $g_\mathrm{f}(V_\mathrm{th})$, which corresponds to a local fast positive feedback, indicates that the neuron is able to fire a spike spontaneously around threshold voltage. A positive $g_\mathrm{s}(V_\mathrm{th})$, which corresponds to a slow negative feedback, indicates that, right after a spike, the neuron will tend to attenuate the excitation and bring back the neuron to rest voltage, while a negative $g_\mathrm{s}(V_\mathrm{th})$ indicates that the neuron will tend to fire other spikes to initiate a burst. In the case of bursting neuron, $g_\mathrm{u}(V_\mathrm{th})$ is always positive and is an indicator of the interburst frequency as well as the duty cycle, i.e., ultraslow negative feedback.

Once the conductance based model had been initialized using *NmodController.jl*, computing DICs or sensitivity matrix can be made using the function `computeDICs()`. As arguments, you only have to specify your model and three time constant functions corresponding to the fast, slow and ultraslow time ranges in which DICs differ to compute these, such that you avoid all the fancy computations by hand.

### Example 1
The next few lines of code show how to compute DICs using a predefined conductance based model `neuron`. Note that you must specify maximum ion channel conductances when initializing the conductance based model in order to compute DICs.

```julia
# First initializing the three timescales
tauX(V, A, B, D, E) = A - B / (1 + exp((V+D) / E))
tauFast(V) = tauX(V, 1.32, 1.26, 120., -25.)
tauSlow(V) = boltz(V, 12.3, -11.8)
tauUltraslow(V) = 100.

using NmodController
gf, gs, gu = computeDICs(neuron, tauFast, tauSlow, tauUltraslow)
```

### Example 2
The next few lines of code show how to compute the sensitivity matrix using a predefined conductance based model `neuron`.

```julia
# First initializing the three timescales
tauX(V, A, B, D, E) = A - B / (1 + exp((V+D) / E))
tauFast(V) = tauX(V, 1.32, 1.26, 120., -25.)
tauSlow(V) = boltz(V, 12.3, -11.8)
tauUltraslow(V) = 100.

using NmodController
S = computeDICs(neuron, tauFast, tauSlow, tauUltraslow, onlyS=true)
```

Note that two other optional arguments can be specified to `computeDICs()`:
1. `tauCa` which is a number capturing the timescale of the intracellular calcium dynamics (often large since calcium is ultraslow);
2. `scaled` which is a boolean automatically set to true that makes DICs or sensitivity matrix being scaled by the leakage conductance of the conductances based model so that these are dimensionless, strongly recommend.

## Computing threshold voltage

As said earlier, threshold voltage $V_\mathrm{th}$ is ubiquitous when dealing with DICs. To compute it, we need DICs as threshold voltage is defined as the positive to negative zero crossing of the global input conductance $g_\mathrm{in}(V) = g_\mathrm{f}(V) + g_\mathrm{s}(V) + g_\mathrm{u}(V)$. Once the DICs had been computed using *NmodController.jl*, computing threshold voltage can be made using the function `computeThresholdVoltage()`.

### Example 3
The next few lines of code show how to compute the threshold voltage using a predefined conductance based model `neuron`.

```julia
# First initializing the three timescales
tauX(V, A, B, D, E) = A - B / (1 + exp((V+D) / E))
tauFast(V) = tauX(V, 1.32, 1.26, 120., -25.)
tauSlow(V) = boltz(V, 12.3, -11.8)
tauUltraslow(V) = 100.

using NmodController
gf, gs, gu = computeDICs(neuron, tauFast, tauSlow, tauUltraslow)
Vth = computeThresholdVoltage(gf, gs, gu)
```