# Initializing Nernst reversalPotential
VNa = 50. # Sodium reversal potential
VK = -80. # Potassium reversal potential
VCa = 80. # Calcium reversal potential
VH = -20. # Reversal potential for the H-current (permeable to both sodium and potassium ions)
Vleak = -50. # Reversal potential of leak channels

# Building Na current
mNa_inf(V) = boltz(V, 25.5, -5.29)
tau_mNa(V) = tauX(V, 1.32, 1.26, 120., -25.)
hNa_inf(V) = boltz(V, 48.9, 5.18)
tau_hNa(V) = (0.67 / (1 + exp((V+62.9)/-10.0))) * (1.5 + 1 / (1+exp((V+34.9)/3.6)))
NaCurrent = initializeCurrent("Na", VNa, numberOfGatings=2, exponents=[3, 1],
    activationSteadyStateGating=mNa_inf, activationTimeConstant=tau_mNa,
    inactivationSteadyStateGating=hNa_inf, inactivationTimeConstant=tau_hNa)

# Building Kd current
mKd_inf(V) = boltz(V, 12.3, -11.8)
tau_mKd(V) = tauX(V, 7.2, 6.4, 28.3, -19.2)
KdCurrent = initializeCurrent("Kd", VK, exponents=4,
    activationSteadyStateGating=mKd_inf, activationTimeConstant=tau_mKd)

# Building KCa current
mKCa_inf(V, Ca) = (Ca / (Ca+3)) * (1 / (1+exp((V+28.3)/-12.6)))
tau_mKCa(V) = tauX(V, 90.3, 75.1, 46., -22.7)
KCaCurrent = initializeCurrent("KCa", VK, exponents=4,
    activationSteadyStateGating=mKCa_inf, activationTimeConstant=tau_mKCa,
    calciumDependency=true)

# Building CaT current
mCaT_inf(V) = boltz(V, 27.1, -7.2)
tau_mCaT(V) = tauX(V, 21.7, 21.3, 68.1, -20.5)
hCaT_inf(V) = boltz(V, 32.1, 5.5)
tau_hCaT(V) = tauX(V, 105., 89.8, 55., -16.9)
CaTCurrent = initializeCurrent("CaT", VCa, numberOfGatings=2, exponents=[3, 1],
    activationSteadyStateGating=mCaT_inf, activationTimeConstant=tau_mCaT,
    inactivationSteadyStateGating=hCaT_inf, inactivationTimeConstant=tau_hCaT)

# Building CaS current
mCaS_inf(V) = boltz(V, 33., -8.1)
tau_mCaS(V) = 1.4 + (7 / ((exp((V+27)/10)) + (exp((V+70)/-13))))
hCaS_inf(V) = boltz(V, 60., 6.2)
tau_hCaS(V) = 60 + (150 / ((exp((V+55)/9)) + (exp((V+65)/-16))))
CaSCurrent = initializeCurrent("CaS", VCa, numberOfGatings=2, exponents=[3, 1],
    activationSteadyStateGating=mCaS_inf, activationTimeConstant=tau_mCaS,
    inactivationSteadyStateGating=hCaS_inf, inactivationTimeConstant=tau_hCaS)

# Building A current
mA_inf(V) = boltz(V, 27.2, -8.7)
tau_mA(V) = tauX(V, 11.6, 10.4, 32.9, -15.2)
hA_inf(V) = boltz(V, 56.9, 4.9)
tau_hA(V) = tauX(V, 38.6, 29.2, 38.9, -26.5)
ACurrent = initializeCurrent("A", VK, numberOfGatings=2, exponents=[3, 1],
    activationSteadyStateGating=mA_inf, activationTimeConstant=tau_mA,
    inactivationSteadyStateGating=hA_inf, inactivationTimeConstant=tau_hA)

# Building H current
mH_inf(V) = boltz(V, 70., 6.)
tau_mH(V) = tauX(V, 272., -1499., 42.2, -8.73)
HCurrent = initializeCurrent("H", VH, exponents=1,
    activationSteadyStateGating=mH_inf, activationTimeConstant=tau_mH)

# Building calcium dynamics
CaDyn = initializeCalciumDynamics(["CaT", "CaS"], [-0.94, -0.94], 0.05, 20)

# Building a more complex model with calcium
ionCurrents = [NaCurrent, CaTCurrent, CaSCurrent, ACurrent, KCaCurrent, KdCurrent, HCurrent]
gvec = [800., 3., 3., 80., 60., 90., 0.1]
STG = initializeNeuronModel(ionCurrents, calciumDynamics=CaDyn, leakageConductance=0.01, reversaleLeakagePotential=Vleak, maximumConductances=gvec)

# Defining some timescales
tauFast = tau_mNa
tauSlow = tau_mKd
tauUltraslow = tau_mH

# Computing S
S = computeDICs(STG, tauFast, tauSlow, tauUltraslow, tauCa=500000., onlyS=true)

# Computing DICs
gf, gs, gu = computeDICs(STG, tauFast, tauSlow, tauUltraslow, tauCa=500000.)

# Returning the test value
tp = 1.25
return (isa(S(tp), Matrix{Float64}) && size(S(tp)) == (3, 7) && isa(gf(tp), Float64) && isa(gs(tp), Float64) && isa(gu(tp), Float64))
