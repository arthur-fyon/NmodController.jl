### STG model
# Building Na current
STG_NaCurrent = initializeCurrent("Na", STG_VNa, numberOfGatings=2, exponents=[3, 1],
    activationSteadyStateGating=STG_mNa_inf, activationTimeConstant=STG_tau_mNa,
    inactivationSteadyStateGating=STG_hNa_inf, inactivationTimeConstant=STG_tau_hNa)

# Building Kd current
STG_KdCurrent = initializeCurrent("Kd", STG_VK, exponents=4,
    activationSteadyStateGating=STG_mKd_inf, activationTimeConstant=STG_tau_mKd)

# Building KCa current
STG_KCaCurrent = initializeCurrent("KCa", STG_VK, exponents=4,
    activationSteadyStateGating=STG_mKCa_inf, activationTimeConstant=STG_tau_mKCa,
    calciumDependency=true)
STG_KCaCurrent2 = initializeCurrent("KCa", STG_VK, numberOfGatings=2, exponents=[3, 1],
    activationSteadyStateGating=STG_mKCa_inf, activationTimeConstant=STG_tau_mKCa,
    inactivationSteadyStateGating=STG_mKCa_inf, inactivationTimeConstant=STG_tau_mKCa,
    calciumDependency=true)

# Building CaT current
STG_CaTCurrent = initializeCurrent("CaT", STG_VCa, numberOfGatings=2, exponents=[3, 1],
    activationSteadyStateGating=STG_mCaT_inf, activationTimeConstant=STG_tau_mCaT,
    inactivationSteadyStateGating=STG_hCaT_inf, inactivationTimeConstant=STG_tau_hCaT)

# Building CaS current
STG_CaSCurrent = initializeCurrent("CaS", STG_VCa, numberOfGatings=2, exponents=[3, 1],
    activationSteadyStateGating=STG_mCaS_inf, activationTimeConstant=STG_tau_mCaS,
    inactivationSteadyStateGating=STG_hCaS_inf, inactivationTimeConstant=STG_tau_hCaS)

# Building A current
STG_ACurrent = initializeCurrent("A", STG_VK, numberOfGatings=2, exponents=[3, 1],
    activationSteadyStateGating=STG_mA_inf, activationTimeConstant=STG_tau_mA,
    inactivationSteadyStateGating=STG_hA_inf, inactivationTimeConstant=STG_tau_hA)

# Building H current
STG_HCurrent = initializeCurrent("H", STG_VH, exponents=1,
    activationSteadyStateGating=STG_mH_inf, activationTimeConstant=STG_tau_mH)

# Building calcium dynamics
CaDyn = initializeCalciumDynamics(["CaT", "CaS"], [-0.94, -0.94], 0.05, 20)
SimpleCaDyn = initializeCalciumDynamics("CaT", -0.94, 0.05, 20)

# Building a more complex model with calcium
STG_ionCurrents = [STG_NaCurrent, STG_CaTCurrent, STG_CaSCurrent, STG_ACurrent, STG_KCaCurrent, STG_KdCurrent, STG_HCurrent]
STG_ionCurrents2 = [STG_NaCurrent, STG_CaTCurrent, STG_CaSCurrent, STG_ACurrent, STG_KCaCurrent2, STG_KdCurrent, STG_HCurrent]
STG_gvec = [800., 3., 3., 80., 60., 90., 0.1]
STG = initializeNeuronModel(STG_ionCurrents, C=0.1, calciumDynamics=CaDyn, leakageConductance=0.01, reversaleLeakagePotential=STG_Vleak, maximumConductances=STG_gvec)
SimpleSTG = initializeNeuronModel(STG_ionCurrents2, calciumDynamics=SimpleCaDyn, maximumConductances=STG_gvec)

### DA model
# Building Na current
DA_NaCurrent = initializeCurrent("Na", DA_VNa, numberOfGatings=2, exponents=[3, 1],
    activationSteadyStateGating=DA_mNa_inf, activationTimeConstant=DA_tau_mNa,
    inactivationSteadyStateGating=DA_hNa_inf, inactivationTimeConstant=DA_tau_hNa)

# Building Kd current
DA_KdCurrent = initializeCurrent("Kd", DA_VK, exponents=3,
    activationSteadyStateGating=DA_mKd_inf, activationTimeConstant=DA_tau_mKd)

# Building CaL current
DA_CaLCurrent = initializeCurrent("CaL", DA_VCa, exponents=2,
    activationSteadyStateGating=DA_mCaL_inf, activationTimeConstant=DA_tau_mCaL)

# Building CaN current
DA_CaSCurrent = initializeCurrent("CaN", DA_VCa, exponents=1,
    activationSteadyStateGating=DA_mCaN_inf, activationTimeConstant=DA_tau_mCaN)

# Building ERG current
DA_ERGCurrent = initializeCurrent("ERG", DA_VK, exponents=1,
    activationSteadyStateGating=DA_o_inf, activationTimeConstant=DA_tau_o)

# Building NMDA current
DA_NMDACurrent = initializeCurrent("NMDA", DA_VNMDA, exponents=1,
    activationSteadyStateGating=DA_NMDA_inf, activationTimeConstant=DA_tau_NMDA,
    MgDependency=true)
DA_NMDACurrent2 = initializeCurrent("NMDA", DA_VNMDA, numberOfGatings=2, exponents=[1, 1],
    activationSteadyStateGating=DA_NMDA_inf, activationTimeConstant=DA_tau_NMDA,
    inactivationSteadyStateGating=DA_NMDA_inf, inactivationTimeConstant=DA_tau_NMDA,
    MgDependency=true)

# Building a more complex model with calcium
DA_ionCurrents = [DA_NaCurrent, DA_KdCurrent, DA_CaLCurrent, DA_CaSCurrent, DA_ERGCurrent, DA_NMDACurrent]
DA_ionCurrents2 = [DA_NaCurrent, DA_KdCurrent, DA_CaLCurrent, DA_CaSCurrent, DA_ERGCurrent, DA_NMDACurrent2]
DA = initializeNeuronModel(DA_ionCurrents, C=1.)
DA2 = initializeNeuronModel(DA_ionCurrents2, C=1.)

# Writing uncontrolled and controlled simulation files
writeUncontrolledODEs(STG)
writeControlledODEs(STG, ["CaS", "A"], ["s", "u"])

# Computing return flag
flag = false
if isfile("CBModelODEs.jl") && isfile("ControlledCBModelODEs.jl")
    rm("CBModelODEs.jl")
    rm("ControlledCBModelODEs.jl")
    flag = true
end

# Writing uncontrolled and controlled simulation files
writeUncontrolledODEs(DA)
writeControlledODEs(DA, ["CaN", "CaL"], ["s"])

# Computing return flag
if isfile("CBModelODEs.jl") && isfile("ControlledCBModelODEs.jl")
    rm("CBModelODEs.jl")
    rm("ControlledCBModelODEs.jl")
    flag = true
else
    flag = false
end

return flag
