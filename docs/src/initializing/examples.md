# Example of existing models in the litterature

Through this documentation, two conductance based model coming from the litterature wil be used with *NmodController.jl*: a stomatogastric (STG) neuron model (Liu et al., 1998 "A model neuron with activity-dependent conductances regulated by multiple calcium sensors") and an adapted dopaminergic (DA) neuron model (Drion et al., 2011 "How modeling can reconcile apparently discrepant experimental results: the case of pacemaking in dopaminergic neurons").

## Initializing the STG model
The STG model is composed of 7 voltage gated ionic currents which one is calcium dependent:
1. Transient sodium current $I_\mathrm{Na}$ (2 gating variables);
2. T-type calcium current $I_\mathrm{CaT}$ (2 gating variables);
3. Slow calcium current $I_\mathrm{CaS}$ (2 gating variables);
4. A-type potassium current $I_\mathrm{A}$ (2 gating variables);
5. Calcium controlled potassium current $I_\mathrm{KCa}$ (1 gating variable caclium dpendent);
6. Delayed rectified potassium current $I_\mathrm{Kd}$ (1 gating variable);
7. H current $I_\mathrm{H}$ (1 gating variable).

The voltage equation writes

$C \dot V = - \bar{g}_\mathrm{Na}m^3_\mathrm{Na}h_\mathrm{Na}(V-E_\mathrm{Na}) - \bar{g}_\mathrm{CaT}m^3_\mathrm{CaT}h_\mathrm{CaT}(V-E_\mathrm{Ca}) - \bar{g}_\mathrm{CaS}m^3_\mathrm{CaS}h_\mathrm{CaS}(V-E_\mathrm{Ca}) - \bar{g}_\mathrm{A}m^3_\mathrm{A}h_\mathrm{A}(V-E_\mathrm{K}) - \bar{g}_\mathrm{KCa}m^4_\mathrm{KCa}(V-E_\mathrm{K}) - \bar{g}_\mathrm{Kd}m^4_\mathrm{Kd}(V-E_\mathrm{K}) - \bar{g}_\mathrm{H}m_\mathrm{H}(V-E_\mathrm{H}) - g_\mathrm{leak}(V-E_\mathrm{leak}) + I_{ext}(t)$

and the intracellular calcium dynamic writes

$\tau_{Ca} \cdot \dot{[Ca]} = -0.94\cdot  I_{CaT} -0.94\cdot  I_{CaS} - [Ca] + Ca_\mathrm{eq}$

with $\tau_{Ca} = 20$ and $Ca_\mathrm{eq} = 0.05$.

The next few lines of code show how to initialize such model using *NmodController.jl*.

```julia
### Kinetics and parameters for the STG model
# STG gating Functions
STG_boltz(V, A, B) = 1 / (1 + exp((V+A) / B))
tauX(V, A, B, D, E) = A - B / (1 + exp((V+D) / E))

# Initializing Nernst reversal potentials
STG_ENa = 50. # Sodium reversal potential
STG_EK = -80. # Potassium reversal potential
STG_ECa = 80. # Calcium reversal potential
STG_EH = -20. # Reversal potential for the H-current (permeable to both sodium and potassium ions)
STG_Eleak = -50. # Reversal potential of leak channels

# Na current
STG_mNa_inf(V) = STG_boltz(V, 25.5, -5.29)
STG_tau_mNa(V) = tauX(V, 1.32, 1.26, 120., -25.)
STG_hNa_inf(V) = STG_boltz(V, 48.9, 5.18)
STG_tau_hNa(V) = (0.67 / (1 + exp((V+62.9)/-10.0))) * (1.5 + 1 / (1+exp((V+34.9)/3.6)))

# Kd current
STG_mKd_inf(V) = STG_boltz(V, 12.3, -11.8)
STG_tau_mKd(V) = tauX(V, 7.2, 6.4, 28.3, -19.2)

# KCa current
STG_mKCa_inf(V, Ca) = (Ca / (Ca+3)) * (1 / (1+exp((V+28.3)/-12.6)))
STG_tau_mKCa(V) = tauX(V, 90.3, 75.1, 46., -22.7)

# CaT current
STG_mCaT_inf(V) = STG_boltz(V, 27.1, -7.2)
STG_tau_mCaT(V) = tauX(V, 21.7, 21.3, 68.1, -20.5)
STG_hCaT_inf(V) = STG_boltz(V, 32.1, 5.5)
STG_tau_hCaT(V) = tauX(V, 105., 89.8, 55., -16.9)

# CaS current
STG_mCaS_inf(V) = STG_boltz(V, 33., -8.1)
STG_tau_mCaS(V) = 1.4 + (7 / ((exp((V+27)/10)) + (exp((V+70)/-13))))
STG_hCaS_inf(V) = STG_boltz(V, 60., 6.2)
STG_tau_hCaS(V) = 60 + (150 / ((exp((V+55)/9)) + (exp((V+65)/-16))))

# A current
STG_mA_inf(V) = STG_boltz(V, 27.2, -8.7)
STG_tau_mA(V) = tauX(V, 11.6, 10.4, 32.9, -15.2)
STG_hA_inf(V) = STG_boltz(V, 56.9, 4.9)
STG_tau_hA(V) = tauX(V, 38.6, 29.2, 38.9, -26.5)

# H current
STG_mH_inf(V) = STG_boltz(V, 70., 6.)
STG_tau_mH(V) = tauX(V, 272., -1499., 42.2, -8.73)

using NmodController

# Building Na current
STG_NaCurrent = initializeCurrent("Na", STG_ENa, numberOfGatings=2, exponents=[3, 1],
    activationSteadyStateGating=STG_mNa_inf, activationTimeConstant=STG_tau_mNa,
    inactivationSteadyStateGating=STG_hNa_inf, inactivationTimeConstant=STG_tau_hNa)

# Building Kd current
STG_KdCurrent = initializeCurrent("Kd", STG_EK, exponents=4,
    activationSteadyStateGating=STG_mKd_inf, activationTimeConstant=STG_tau_mKd)

# Building KCa current
STG_KCaCurrent = initializeCurrent("KCa", STG_EK, exponents=4,
    activationSteadyStateGating=STG_mKCa_inf, activationTimeConstant=STG_tau_mKCa,
    calciumDependency=true)

# Building CaT current
STG_CaTCurrent = initializeCurrent("CaT", STG_ECa, numberOfGatings=2, exponents=[3, 1],
    activationSteadyStateGating=STG_mCaT_inf, activationTimeConstant=STG_tau_mCaT,
    inactivationSteadyStateGating=STG_hCaT_inf, inactivationTimeConstant=STG_tau_hCaT)

# Building CaS current
STG_CaSCurrent = initializeCurrent("CaS", STG_ECa, numberOfGatings=2, exponents=[3, 1],
    activationSteadyStateGating=STG_mCaS_inf, activationTimeConstant=STG_tau_mCaS,
    inactivationSteadyStateGating=STG_hCaS_inf, inactivationTimeConstant=STG_tau_hCaS)

# Building A current
STG_ACurrent = initializeCurrent("A", STG_EK, numberOfGatings=2, exponents=[3, 1],
    activationSteadyStateGating=STG_mA_inf, activationTimeConstant=STG_tau_mA,
    inactivationSteadyStateGating=STG_hA_inf, inactivationTimeConstant=STG_tau_hA)

# Building H current
STG_HCurrent = initializeCurrent("H", STG_EH, exponents=1,
    activationSteadyStateGating=STG_mH_inf, activationTimeConstant=STG_tau_mH)

# Building calcium dynamics
CaDyn = initializeCalciumDynamics(["CaT", "CaS"], [-0.94, -0.94], 0.05, 20)

# Wrapping all currents in a vector
STG_ionCurrents = [STG_NaCurrent, STG_CaTCurrent, STG_CaSCurrent, STG_ACurrent, STG_KCaCurrent, STG_KdCurrent, STG_HCurrent]
STG_gvec = [800., 3., 3., 80., 60., 90., 0.1]

# Initializing the STG model
STG = initializeNeuronModel(STG_ionCurrents, C=0.1, calciumDynamics=CaDyn, leakageConductance=0.01, reversaleLeakagePotential=STG_Eleak, maximumConductances=STG_gvec)
```