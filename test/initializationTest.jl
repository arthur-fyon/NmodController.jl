# Building Na current
mNa_inf(V) = boltz(V, 25.5, -5.29)
tau_mNa(V) = tauX(V, 1.32, 1.26, 120., -25.)
hNa_inf(V) = boltz(V, 48.9, 5.18)
tau_hNa(V) = (0.67 / (1 + exp((V+62.9)/-10.0))) * (1.5 + 1 / (1+exp((V+34.9)/3.6)))
NaCurrent = initializeCurrent("Na", 20, numberOfGatings=2, exponents=[3, 1],
    activationSteadyStateGating=mNa_inf, activationTimeConstants=tau_mNa,
    inactivationSteadyStateGating=hNa_inf, inactivationTimeConstants=tau_hNa)

# Building Kd current
mKd_inf(V) = boltz(V, 12.3, -11.8)
tau_mKd(V) = tauX(V, 7.2, 6.4, 28.3, -19.2)
KdCurrent = initializeCurrent("Kd", -50, exponents=4,
    activationSteadyStateGating=mKd_inf, activationTimeConstants=tau_mKd)

# Building Hodgkin and Huxley model
ionCurrentsHH = [NaCurrent, KdCurrent]
HHmodel = initializeNeuronModel(ionCurrentsHH, C=0.1, maximumConductances=[200., 50.])

# Building KCa current
mKCa_inf(V, Ca) = (Ca / (Ca+3)) * (1 / (1+exp((V+28.3)/-12.6)))
tau_mKCa(V) = tauX(V, 90.3, 75.1, 46., -22.7)
KCaCurrent = initializeCurrent("KCa", -50, exponents=4,
    activationSteadyStateGating=mKCa_inf, activationTimeConstants=tau_mKCa,
    calciumDependency=true)

# Building CaT current
mCaT_inf(V) = boltz(V, 27.1, -7.2)
tau_mCaT(V) = tauX(V, 21.7, 21.3, 68.1, -20.5)
hCaT_inf(V) = boltz(V, 32.1, 5.5)
tau_hCaT(V) = tauX(V, 105., 89.8, 55., -16.9)
CaTCurrent = initializeCurrent("CaT", -70, numberOfGatings=2, exponents=[3, 1],
    activationSteadyStateGating=mCaT_inf, activationTimeConstants=tau_mCaT,
    inactivationSteadyStateGating=hCaT_inf, inactivationTimeConstants=tau_hCaT)

# Building CaS current
mCaS_inf(V) = boltz(V, 33., -8.1)
tau_mCaS(V) = 1.4 + (7 / ((exp((V+27)/10)) + (exp((V+70)/-13))))
hCaS_inf(V) = boltz(V, 60., 6.2)
tau_hCaS(V) = 60 + (150 / ((exp((V+55)/9)) + (exp((V+65)/-16))))
CaSCurrent = initializeCurrent("CaS", -70, numberOfGatings=2, exponents=[3, 1],
    activationSteadyStateGating=mCaS_inf, activationTimeConstants=tau_mCaS,
    inactivationSteadyStateGating=hCaS_inf, inactivationTimeConstants=tau_hCaS)

# Building calcium dynamics
CaDyn = initializeCalciumDynamics(["CaT", "CaS"], [-0.94, -0.94], 0.05, 20)

# Building a more complex model with calcium
ionCurrents = [NaCurrent, CaTCurrent, KCaCurrent, CaSCurrent, KdCurrent]
bigModel = initializeNeuronModel(ionCurrents, C=0.1, calciumDynamics=CaDyn)

return (isa(CaTCurrent, IonCurrent) && isa(CaDyn, CalciumDynamic) && isa(HHmodel, NeuronCB) && isa(bigModel, NeuronCB))
