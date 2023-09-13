var documenterSearchIndex = {"docs":
[{"location":"computing/DICs/#Computations-on-conductance-based-model","page":"How to compute DICs and sensitivity matrix?","title":"Computations on conductance based model","text":"","category":"section"},{"location":"computing/DICs/#Computing-dynamic-input-conductances-(DICs)","page":"How to compute DICs and sensitivity matrix?","title":"Computing dynamic input conductances (DICs)","text":"","category":"section"},{"location":"computing/DICs/","page":"How to compute DICs and sensitivity matrix?","title":"How to compute DICs and sensitivity matrix?","text":"DICs are three voltage-dependent conductance curves g_mathrmf(V)g_mathrms(V)g_mathrmu(V) that can  be computed as linear functions of the maximal conductance vector barg_mathrmion of the neuron model at each V. By wrapping the three dynamic input conductances in a vector g_mathrmDIC(V), the linear combination writes","category":"page"},{"location":"computing/DICs/","page":"How to compute DICs and sensitivity matrix?","title":"How to compute DICs and sensitivity matrix?","text":"g_mathrmDIC(V) = S(V) cdot barg_mathrmion","category":"page"},{"location":"computing/DICs/","page":"How to compute DICs and sensitivity matrix?","title":"How to compute DICs and sensitivity matrix?","text":"where S(V) is a sensitivity matrix. These mathematical function can be built using Drion et al., 2015 \"Dynamic input conductances shape neuronal spiking\".","category":"page"},{"location":"computing/DICs/","page":"How to compute DICs and sensitivity matrix?","title":"How to compute DICs and sensitivity matrix?","text":"Because of the specific feedback structure of conductance-based models, DICs shape neuronal spiking behavior and the three DICs differ in the timescale at which this shaping happens: fast g_mathrmf(V), slow g_mathrms(V), and ultraslow g_mathrmu(V). Values and signs of the DICs at specific voltages, mainly the threshold voltage V_mathrmth, reliably determine the neuronal firing pattern. For instance, a negative g_mathrmf(V_mathrmth), which corresponds to a local fast positive feedback, indicates that the neuron is able to fire a spike spontaneously around threshold voltage. A positive g_mathrms(V_mathrmth), which corresponds to a slow negative feedback, indicates that, right after a spike, the neuron will tend to attenuate the excitation and bring back the neuron to rest voltage, while a negative g_mathrms(V_mathrmth) indicates that the neuron will tend to fire other spikes to initiate a burst. In the case of bursting neuron, g_mathrmu(V_mathrmth) is always positive and is an indicator of the interburst frequency as well as the duty cycle, i.e., ultraslow negative feedback.","category":"page"},{"location":"computing/DICs/","page":"How to compute DICs and sensitivity matrix?","title":"How to compute DICs and sensitivity matrix?","text":"Once the conductance based model had been initialized using NmodController.jl, computing DICs or sensitivity matrix can be made using the function computeDICs(). As arguments, you only have to specify your model and three time constant functions corresponding to the fast, slow and ultraslow time ranges in which DICs differ to compute these, such that you avoid all the fancy computations by hand.","category":"page"},{"location":"computing/DICs/#Example-1","page":"How to compute DICs and sensitivity matrix?","title":"Example 1","text":"","category":"section"},{"location":"computing/DICs/","page":"How to compute DICs and sensitivity matrix?","title":"How to compute DICs and sensitivity matrix?","text":"The next few lines of code show how to compute DICs using a predefined conductance based model neuron. Note that you must specify maximum ion channel conductances when initializing the conductance based model in order to compute DICs.","category":"page"},{"location":"computing/DICs/","page":"How to compute DICs and sensitivity matrix?","title":"How to compute DICs and sensitivity matrix?","text":"# First initializing the three timescales\ntauX(V, A, B, D, E) = A - B / (1 + exp((V+D) / E))\ntauFast(V) = tauX(V, 1.32, 1.26, 120., -25.)\ntauSlow(V) = boltz(V, 12.3, -11.8)\ntauUltraslow(V) = 100.\n\nusing NmodController\ngf, gs, gu = computeDICs(neuron, tauFast, tauSlow, tauUltraslow)","category":"page"},{"location":"computing/DICs/#Example-2","page":"How to compute DICs and sensitivity matrix?","title":"Example 2","text":"","category":"section"},{"location":"computing/DICs/","page":"How to compute DICs and sensitivity matrix?","title":"How to compute DICs and sensitivity matrix?","text":"The next few lines of code show how to compute the sensitivity matrix using a predefined conductance based model neuron.","category":"page"},{"location":"computing/DICs/","page":"How to compute DICs and sensitivity matrix?","title":"How to compute DICs and sensitivity matrix?","text":"# First initializing the three timescales\ntauX(V, A, B, D, E) = A - B / (1 + exp((V+D) / E))\ntauFast(V) = tauX(V, 1.32, 1.26, 120., -25.)\ntauSlow(V) = boltz(V, 12.3, -11.8)\ntauUltraslow(V) = 100.\n\nusing NmodController\nS = computeDICs(neuron, tauFast, tauSlow, tauUltraslow, onlyS=true)","category":"page"},{"location":"computing/DICs/","page":"How to compute DICs and sensitivity matrix?","title":"How to compute DICs and sensitivity matrix?","text":"Note that two other optional arguments can be specified to computeDICs():","category":"page"},{"location":"computing/DICs/","page":"How to compute DICs and sensitivity matrix?","title":"How to compute DICs and sensitivity matrix?","text":"tauCa which is a number capturing the timescale of the intracellular calcium dynamics (often large since calcium is ultraslow);\nscaled which is a boolean automatically set to true that makes DICs or sensitivity matrix being scaled by the leakage conductance of the conductances based model so that these are dimensionless, strongly recommend.","category":"page"},{"location":"computing/DICs/#Computing-threshold-voltage","page":"How to compute DICs and sensitivity matrix?","title":"Computing threshold voltage","text":"","category":"section"},{"location":"computing/DICs/","page":"How to compute DICs and sensitivity matrix?","title":"How to compute DICs and sensitivity matrix?","text":"As said earlier, threshold voltage V_mathrmth is ubiquitous when dealing with DICs. To compute it, we need DICs as threshold voltage is defined as the positive to negative zero crossing of the global input conductance g_mathrmin(V) = g_mathrmf(V) + g_mathrms(V) + g_mathrmu(V). Once the DICs had been computed using NmodController.jl, computing threshold voltage can be made using the function computeThresholdVoltage().","category":"page"},{"location":"computing/DICs/#Example-3","page":"How to compute DICs and sensitivity matrix?","title":"Example 3","text":"","category":"section"},{"location":"computing/DICs/","page":"How to compute DICs and sensitivity matrix?","title":"How to compute DICs and sensitivity matrix?","text":"The next few lines of code show how to compute the threshold voltage using a predefined conductance based model neuron.","category":"page"},{"location":"computing/DICs/","page":"How to compute DICs and sensitivity matrix?","title":"How to compute DICs and sensitivity matrix?","text":"# First initializing the three timescales\ntauX(V, A, B, D, E) = A - B / (1 + exp((V+D) / E))\ntauFast(V) = tauX(V, 1.32, 1.26, 120., -25.)\ntauSlow(V) = boltz(V, 12.3, -11.8)\ntauUltraslow(V) = 100.\n\nusing NmodController\ngf, gs, gu = computeDICs(neuron, tauFast, tauSlow, tauUltraslow)\nVth = computeThresholdVoltage(gf, gs, gu)","category":"page"},{"location":"initializing/current/#Initializing-an-ionic-current","page":"How to initialize an ionic current?","title":"Initializing an ionic current","text":"","category":"section"},{"location":"initializing/current/","page":"How to initialize an ionic current?","title":"How to initialize an ionic current?","text":"In conductance based models, an ionic current I_mathrmion is described by its voltage dependent conductance g_mathrmion and Nernst reversal potential E_mathrmion","category":"page"},{"location":"initializing/current/","page":"How to initialize an ionic current?","title":"How to initialize an ionic current?","text":"I_mathrmion = g_mathrmion(V) cdot (V - E_mathrmion)","category":"page"},{"location":"initializing/current/","page":"How to initialize an ionic current?","title":"How to initialize an ionic current?","text":"This voltage dependent conductance is defined by its maximum value barg_mathrmion and at most two gating variables, one activation m_mathrmion(V) and one inactivation h_mathrmion(V), that varies between 0 and 1. In the case of a current with two gating variables, this writes","category":"page"},{"location":"initializing/current/","page":"How to initialize an ionic current?","title":"How to initialize an ionic current?","text":"g_mathrmion(V) = barg_mathrmion cdot m^a_mathrmion(V) cdot h^b_mathrmion(V)","category":"page"},{"location":"initializing/current/","page":"How to initialize an ionic current?","title":"How to initialize an ionic current?","text":"Note that the maximum ion channel conductance of the current will not be stocked inside the current data structure, but rather in the conductance based model data structure.","category":"page"},{"location":"initializing/current/","page":"How to initialize an ionic current?","title":"How to initialize an ionic current?","text":"Each gating variable dynamic follow a basic first order ODE where both the time constant tau_m_mathrmion and the converging value m_mathrmion_infty are voltage dependent: ","category":"page"},{"location":"initializing/current/","page":"How to initialize an ionic current?","title":"How to initialize an ionic current?","text":"dotm_mathrmion = fracm_mathrmion_infty(V) - m_mathrmiontau_m_mathrmion(V)","category":"page"},{"location":"initializing/current/","page":"How to initialize an ionic current?","title":"How to initialize an ionic current?","text":"In NmodController.jl, an ionic current can be contained in a IonCurrent type. To help initializing such data structure, calling initializeCurrent() with appropriate arguments is strongly recommended.","category":"page"},{"location":"initializing/current/#Example-1","page":"How to initialize an ionic current?","title":"Example 1","text":"","category":"section"},{"location":"initializing/current/","page":"How to initialize an ionic current?","title":"How to initialize an ionic current?","text":"The next few lines of code show how to initialize a sodium current with two gating variables. This current writes","category":"page"},{"location":"initializing/current/","page":"How to initialize an ionic current?","title":"How to initialize an ionic current?","text":"I_mathrmNa = barg_mathrmNa cdot m^3_mathrmNa cdot h_mathrmNa cdot (V - E_mathrmNa)","category":"page"},{"location":"initializing/current/","page":"How to initialize an ionic current?","title":"How to initialize an ionic current?","text":"# First initializing the converging values and time constants functions\nboltz(V, A, B) = 1 / (1 + exp((V+A) / B))\ntauX(V, A, B, D, E) = A - B / (1 + exp((V+D) / E))\n\nmNa_inf(V) = boltz(V, 25.5, -5.29)\ntau_mNa(V) = tauX(V, 1.32, 1.26, 120., -25.)\nhNa_inf(V) = boltz(V, 48.9, 5.18)\ntau_hNa(V) = (0.67 / (1 + exp((V+62.9)/-10.0))) * (1.5 + 1 / (1+exp((V+34.9)/3.6)))\n\n# Initializing Nernst reversal potential\nENa = 50.\n\nusing NmodController\nNaCurrent = initializeCurrent(\"Na\", ENa, numberOfGatings=2, exponents=[3, 1],\n    activationSteadyStateGating=mNa_inf, activationTimeConstant=tau_mNa,\n    inactivationSteadyStateGating=hNa_inf, inactivationTimeConstant=tau_hNa)","category":"page"},{"location":"initializing/current/#Example-2","page":"How to initialize an ionic current?","title":"Example 2","text":"","category":"section"},{"location":"initializing/current/","page":"How to initialize an ionic current?","title":"How to initialize an ionic current?","text":"Sometimes, the gating variable converging values may depend on the intracellular calcium. The next few lines of code show how to initialize a calcium controlled potassium current with one gating variable. This current writes","category":"page"},{"location":"initializing/current/","page":"How to initialize an ionic current?","title":"How to initialize an ionic current?","text":"I_mathrmKCa = barg_mathrmKCa cdot m^4_mathrmKCa cdot (V - E_mathrmKCa)","category":"page"},{"location":"initializing/current/","page":"How to initialize an ionic current?","title":"How to initialize an ionic current?","text":"Where the dynamic of m_mathrmKCa is described by","category":"page"},{"location":"initializing/current/","page":"How to initialize an ionic current?","title":"How to initialize an ionic current?","text":"dotm_mathrmKCa = fracm_mathrmKCa_infty(V Ca) - m_mathrmKCatau_m_mathrmKCa(V)","category":"page"},{"location":"initializing/current/","page":"How to initialize an ionic current?","title":"How to initialize an ionic current?","text":"As you can see, m_mathrmKCa_infty depends on both the voltage and the calcium, this must be notified to initializeCurrent().","category":"page"},{"location":"initializing/current/","page":"How to initialize an ionic current?","title":"How to initialize an ionic current?","text":"# First initializing the converging value and time constant functions\nboltz(V, A, B) = 1 / (1 + exp((V+A) / B))\ntauX(V, A, B, D, E) = A - B / (1 + exp((V+D) / E))\n\nmKCa_inf(V, Ca) = (Ca / (Ca+3)) * (1 / (1+exp((V+28.3)/-12.6)))\ntau_mKCa(V) = tauX(V, 90.3, 75.1, 46., -22.7)\n\n# Initializing Nernst reversal potential\nEK = -80.\n\nusing NmodController\nKCaCurrent = initializeCurrent(\"KCa\", EK, exponents=4,\n    activationSteadyStateGating=mKCa_inf, activationTimeConstant=tau_mKCa,\n    calciumDependency=true)","category":"page"},{"location":"initializing/current/#Example-3","page":"How to initialize an ionic current?","title":"Example 3","text":"","category":"section"},{"location":"initializing/current/","page":"How to initialize an ionic current?","title":"How to initialize an ionic current?","text":"Another specific case is where the time constant may not depend on the voltage. In such case, just provide the time constant as a Float64 or a Int64 in argument of initializeCurrent(). Moreover, the converging value of the gating variable may be magnesium dependent. In such case, do as the following. The next few lines of code show how to initialize an instantaneous magnesium dependent NMDA current with one gating variable. This current writes","category":"page"},{"location":"initializing/current/","page":"How to initialize an ionic current?","title":"How to initialize an ionic current?","text":"I_mathrmNMDA = barg_mathrmNMDA cdot m_mathrmNMDA cdot (V - E_mathrmNMDA)","category":"page"},{"location":"initializing/current/","page":"How to initialize an ionic current?","title":"How to initialize an ionic current?","text":"Where the dynamic of m_mathrmNMDA is described by","category":"page"},{"location":"initializing/current/","page":"How to initialize an ionic current?","title":"How to initialize an ionic current?","text":"m_mathrmNMDA = m_mathrmNMDA_infty(V Mg)","category":"page"},{"location":"initializing/current/","page":"How to initialize an ionic current?","title":"How to initialize an ionic current?","text":"As you can see, m_mathrmNMDA_infty depends on both the voltage and the magnesium, this must be notified to initializeCurrent()","category":"page"},{"location":"initializing/current/","page":"How to initialize an ionic current?","title":"How to initialize an ionic current?","text":"# First initializing the converging value and time constant functions\nNMDA_inf(V, Mg) = 1 / (1 + Mg*exp(-0.08*V)/10.)\ntau_NMDA(V) = 1e-10 # Never equal to zero, this might be fixed later on\n\n# Initializing Nernst reversal potential\nENMDA = 0.\n\nusing NmodController\nNMDACurrent = initializeCurrent(\"NMDA\", ENMDA, exponents=1,\n    activationSteadyStateGating=NMDA_inf, activationTimeConstant=tau_NMDA,\n    MgDependency=true)","category":"page"},{"location":"computing/examples/#WIP","page":"Examples of existing models","title":"WIP","text":"","category":"section"},{"location":"simulating/examples/#WIP","page":"Examples of existing models","title":"WIP","text":"","category":"section"},{"location":"simulating/files/#WIP","page":"How to write ODE functions file?","title":"WIP","text":"","category":"section"},{"location":"initializing/examples/#Example-of-existing-models-in-the-litterature","page":"Examples of existing models","title":"Example of existing models in the litterature","text":"","category":"section"},{"location":"initializing/examples/","page":"Examples of existing models","title":"Examples of existing models","text":"Through this documentation, two conductance based model coming from the litterature wil be used with NmodController.jl: a stomatogastric (STG) neuron model (Liu et al., 1998 \"A model neuron with activity-dependent conductances regulated by multiple calcium sensors\") and an adapted dopaminergic (DA) neuron model (Drion et al., 2011 \"How modeling can reconcile apparently discrepant experimental results: the case of pacemaking in dopaminergic neurons\").","category":"page"},{"location":"initializing/examples/#Initializing-the-STG-model","page":"Examples of existing models","title":"Initializing the STG model","text":"","category":"section"},{"location":"initializing/examples/","page":"Examples of existing models","title":"Examples of existing models","text":"The STG model is composed of 7 voltage gated ionic currents which one is calcium dependent:","category":"page"},{"location":"initializing/examples/","page":"Examples of existing models","title":"Examples of existing models","text":"Transient sodium current I_mathrmNa (2 gating variables);\nT-type calcium current I_mathrmCaT (2 gating variables);\nSlow calcium current I_mathrmCaS (2 gating variables);\nA-type potassium current I_mathrmA (2 gating variables);\nCalcium controlled potassium current I_mathrmKCa (1 gating variable calcium dpendent);\nDelayed rectified potassium current I_mathrmKd (1 gating variable);\nH current I_mathrmH (1 gating variable).","category":"page"},{"location":"initializing/examples/","page":"Examples of existing models","title":"Examples of existing models","text":"The voltage equation writes","category":"page"},{"location":"initializing/examples/","page":"Examples of existing models","title":"Examples of existing models","text":"C dot V = - barg_mathrmNam^3_mathrmNah_mathrmNa(V-E_mathrmNa) - barg_mathrmCaTm^3_mathrmCaTh_mathrmCaT(V-E_mathrmCa) - barg_mathrmCaSm^3_mathrmCaSh_mathrmCaS(V-E_mathrmCa) - barg_mathrmAm^3_mathrmAh_mathrmA(V-E_mathrmK) - barg_mathrmKCam^4_mathrmKCa(V-E_mathrmK) - barg_mathrmKdm^4_mathrmKd(V-E_mathrmK) - barg_mathrmHm_mathrmH(V-E_mathrmH) - g_mathrmleak(V-E_mathrmleak) + I_ext","category":"page"},{"location":"initializing/examples/","page":"Examples of existing models","title":"Examples of existing models","text":"and the intracellular calcium dynamic writes","category":"page"},{"location":"initializing/examples/","page":"Examples of existing models","title":"Examples of existing models","text":"tau_Ca cdot dotCa = -094cdot  I_CaT -094cdot  I_CaS - Ca + Ca_mathrmeq","category":"page"},{"location":"initializing/examples/","page":"Examples of existing models","title":"Examples of existing models","text":"with tau_Ca = 20 and Ca_mathrmeq = 005.","category":"page"},{"location":"initializing/examples/","page":"Examples of existing models","title":"Examples of existing models","text":"The next few lines of code show how to initialize such model using NmodController.jl.","category":"page"},{"location":"initializing/examples/","page":"Examples of existing models","title":"Examples of existing models","text":"### Kinetics and parameters for the STG model\n# STG gating Functions\nSTG_boltz(V, A, B) = 1 / (1 + exp((V+A) / B))\ntauX(V, A, B, D, E) = A - B / (1 + exp((V+D) / E))\n\n# Initializing Nernst reversal potentials\nSTG_ENa = 50. # Sodium reversal potential\nSTG_EK = -80. # Potassium reversal potential\nSTG_ECa = 80. # Calcium reversal potential\nSTG_EH = -20. # Reversal potential for the H-current (permeable to both sodium and potassium ions)\nSTG_Eleak = -50. # Reversal potential of leak channels\n\n# Na current\nSTG_mNa_inf(V) = STG_boltz(V, 25.5, -5.29)\nSTG_tau_mNa(V) = tauX(V, 1.32, 1.26, 120., -25.)\nSTG_hNa_inf(V) = STG_boltz(V, 48.9, 5.18)\nSTG_tau_hNa(V) = (0.67 / (1 + exp((V+62.9)/-10.0))) * (1.5 + 1 / (1+exp((V+34.9)/3.6)))\n\n# Kd current\nSTG_mKd_inf(V) = STG_boltz(V, 12.3, -11.8)\nSTG_tau_mKd(V) = tauX(V, 7.2, 6.4, 28.3, -19.2)\n\n# KCa current\nSTG_mKCa_inf(V, Ca) = (Ca / (Ca+3)) * (1 / (1+exp((V+28.3)/-12.6)))\nSTG_tau_mKCa(V) = tauX(V, 90.3, 75.1, 46., -22.7)\n\n# CaT current\nSTG_mCaT_inf(V) = STG_boltz(V, 27.1, -7.2)\nSTG_tau_mCaT(V) = tauX(V, 21.7, 21.3, 68.1, -20.5)\nSTG_hCaT_inf(V) = STG_boltz(V, 32.1, 5.5)\nSTG_tau_hCaT(V) = tauX(V, 105., 89.8, 55., -16.9)\n\n# CaS current\nSTG_mCaS_inf(V) = STG_boltz(V, 33., -8.1)\nSTG_tau_mCaS(V) = 1.4 + (7 / ((exp((V+27)/10)) + (exp((V+70)/-13))))\nSTG_hCaS_inf(V) = STG_boltz(V, 60., 6.2)\nSTG_tau_hCaS(V) = 60 + (150 / ((exp((V+55)/9)) + (exp((V+65)/-16))))\n\n# A current\nSTG_mA_inf(V) = STG_boltz(V, 27.2, -8.7)\nSTG_tau_mA(V) = tauX(V, 11.6, 10.4, 32.9, -15.2)\nSTG_hA_inf(V) = STG_boltz(V, 56.9, 4.9)\nSTG_tau_hA(V) = tauX(V, 38.6, 29.2, 38.9, -26.5)\n\n# H current\nSTG_mH_inf(V) = STG_boltz(V, 70., 6.)\nSTG_tau_mH(V) = tauX(V, 272., -1499., 42.2, -8.73)\n\nusing NmodController\n\n# Building Na current\nSTG_NaCurrent = initializeCurrent(\"Na\", STG_ENa, numberOfGatings=2, exponents=[3, 1],\n    activationSteadyStateGating=STG_mNa_inf, activationTimeConstant=STG_tau_mNa,\n    inactivationSteadyStateGating=STG_hNa_inf, inactivationTimeConstant=STG_tau_hNa)\n\n# Building Kd current\nSTG_KdCurrent = initializeCurrent(\"Kd\", STG_EK, exponents=4,\n    activationSteadyStateGating=STG_mKd_inf, activationTimeConstant=STG_tau_mKd)\n\n# Building KCa current\nSTG_KCaCurrent = initializeCurrent(\"KCa\", STG_EK, exponents=4,\n    activationSteadyStateGating=STG_mKCa_inf, activationTimeConstant=STG_tau_mKCa,\n    calciumDependency=true)\n\n# Building CaT current\nSTG_CaTCurrent = initializeCurrent(\"CaT\", STG_ECa, numberOfGatings=2, exponents=[3, 1],\n    activationSteadyStateGating=STG_mCaT_inf, activationTimeConstant=STG_tau_mCaT,\n    inactivationSteadyStateGating=STG_hCaT_inf, inactivationTimeConstant=STG_tau_hCaT)\n\n# Building CaS current\nSTG_CaSCurrent = initializeCurrent(\"CaS\", STG_ECa, numberOfGatings=2, exponents=[3, 1],\n    activationSteadyStateGating=STG_mCaS_inf, activationTimeConstant=STG_tau_mCaS,\n    inactivationSteadyStateGating=STG_hCaS_inf, inactivationTimeConstant=STG_tau_hCaS)\n\n# Building A current\nSTG_ACurrent = initializeCurrent(\"A\", STG_EK, numberOfGatings=2, exponents=[3, 1],\n    activationSteadyStateGating=STG_mA_inf, activationTimeConstant=STG_tau_mA,\n    inactivationSteadyStateGating=STG_hA_inf, inactivationTimeConstant=STG_tau_hA)\n\n# Building H current\nSTG_HCurrent = initializeCurrent(\"H\", STG_EH, exponents=1,\n    activationSteadyStateGating=STG_mH_inf, activationTimeConstant=STG_tau_mH)\n\n# Building calcium dynamics\nCaDyn = initializeCalciumDynamics([\"CaT\", \"CaS\"], [-0.94, -0.94], 0.05, 20)\n\n# Wrapping all currents in a vector\nSTG_ionCurrents = [STG_NaCurrent, STG_CaTCurrent, STG_CaSCurrent, STG_ACurrent, STG_KCaCurrent, STG_KdCurrent, STG_HCurrent]\nSTG_gvec = [800., 3., 3., 80., 60., 90., 0.1]\n\n# Initializing the STG model\nSTG = initializeNeuronModel(STG_ionCurrents, C=0.1, calciumDynamics=CaDyn, leakageConductance=0.01, reversaleLeakagePotential=STG_Eleak, maximumConductances=STG_gvec)","category":"page"},{"location":"initializing/examples/#Initializing-the-DA-model","page":"Examples of existing models","title":"Initializing the DA model","text":"","category":"section"},{"location":"initializing/examples/","page":"Examples of existing models","title":"Examples of existing models","text":"The STG model is composed of 6 voltage gated ionic currents which one is magnesium dependent:","category":"page"},{"location":"initializing/examples/","page":"Examples of existing models","title":"Examples of existing models","text":"Transient sodium current I_mathrmNa (2 gating variables);\nDelayed rectified potassium current I_mathrmKd (1 gating variable);\nL-type calcium current I_mathrmCaL (1 gating variable);\nN-type calcium current I_mathrmCaN (1 gating variable);\nERG potassium current I_mathrmERG (1 gating variable);\nNMDA current I_mathrmNMDA (1 gating variable magnesium dpendent).","category":"page"},{"location":"initializing/examples/","page":"Examples of existing models","title":"Examples of existing models","text":"The voltage equation writes","category":"page"},{"location":"initializing/examples/","page":"Examples of existing models","title":"Examples of existing models","text":"C dot V = - barg_mathrmNam^3_mathrmNah_mathrmNa(V-E_mathrmNa) - barg_mathrmKdm^3_mathrmKd(V-E_mathrmK) - barg_mathrmCaLm^2_mathrmCaL(V-E_mathrmCa) - barg_mathrmCaNm_mathrmCaN(V-E_mathrmCa) - barg_mathrmERGm_mathrmERG(V-E_mathrmK) - barg_mathrmNMDAm_mathrmNMDA(V-E_mathrmNMDA) - g_mathrmleak(V-E_mathrmleak) + I_ext","category":"page"},{"location":"initializing/examples/","page":"Examples of existing models","title":"Examples of existing models","text":"The next few lines of code show how to initialize such model using NmodController.jl.","category":"page"},{"location":"initializing/examples/","page":"Examples of existing models","title":"Examples of existing models","text":"### Kinetics and parameters for the DA model\n# DA gating Functions\nDA_boltz(V, A, B) = 1 / (1 + exp(-(V-A) / B))\n\n# Initializing Nernst reversalPotential\nDA_ENa = 60. # Sodium reversal potential\nDA_EK = -85. # Potassium reversal potential\nDA_ECa = 60. # Calcium reversal potential\nDA_ENMDA = 0. # NMDA reversal potential\nDA_Eleak = -50. # Reversal potential of leak channels\n\n# Na current\nDA_mNa_inf(V) = DA_boltz(V, -30.0907, 9.7264)\nDA_tau_mNa(V) = 0.01 + 1.0 / ((-(15.6504 + 0.4043*V)/(exp(-19.565 -0.5052*V)-1.0)) + 3.0212*exp(-7.4630e-3*V))\nDA_hNa_inf(V) = DA_boltz(V, -54.0289, -10.7665)\nDA_tau_hNa(V) = 0.4 + 1.0 / ((5.0754e-4*exp(-6.3213e-2*V)) + 9.7529*exp(0.13442*V))\n\n# Kd current\nDA_mKd_inf(V) = DA_boltz(V, -25., 12.)\nDA_tau_mKd(V) = tauX(V, 20., 18., 38., -10.)\n\n# CaL current\nDA_mCaL_inf(V) = DA_boltz(V, -50., 2.)\nDA_tau_mCaL(V) = tauX(V, 30., 28., 45., -3.)\n\n# CaN current\nDA_mCaN_inf(V) = DA_boltz(V, -30., 7.)\nDA_tau_mCaN(V) = tauX(V, 30., 25., 55., -6.)\n\n# ERG current (this current is a bit special and its gating dynamic is modified to fit\n# in the generalized form described earlier)\na0ERG(V) = 0.0036 * exp(0.0759*V)\nb0ERG(V) = 1.2523e-5 * exp(-0.0671*V)\naiERG(V) = 0.1 * exp(0.1189*V)\nbiERG(V) = 0.003 * exp(-0.0733*V)\nDA_o_inf(V) = a0ERG(V)*biERG(V) / (a0ERG(V)*(aiERG(V)+biERG(V)) + b0ERG(V)*biERG(V))\n# Due to this generalization, we don't have any time constant. However, this current \n# is ultraslow anyway, so putting an extremely large time constant is the good way to do\nDA_tau_o(V) = 100000. \n\n# NMDA current\nDA_NMDA_inf(V, Mg) = 1 / (1 + Mg*exp(-0.08*V)/10.)\n# This gating variable is instantaneous, however, putting a null time constant is not supported yet\nDA_tau_NMDA(V) = 1e-10\n\nusing NmodController\n\n# Building Na current\nDA_NaCurrent = initializeCurrent(\"Na\", DA_ENa, numberOfGatings=2, exponents=[3, 1],\n    activationSteadyStateGating=DA_mNa_inf, activationTimeConstant=DA_tau_mNa,\n    inactivationSteadyStateGating=DA_hNa_inf, inactivationTimeConstant=DA_tau_hNa)\n\n# Building Kd current\nDA_KdCurrent = initializeCurrent(\"Kd\", DA_EK, exponents=3,\n    activationSteadyStateGating=DA_mKd_inf, activationTimeConstant=DA_tau_mKd)\n\n# Building CaL current\nDA_CaLCurrent = initializeCurrent(\"CaL\", DA_ECa, exponents=2,\n    activationSteadyStateGating=DA_mCaL_inf, activationTimeConstant=DA_tau_mCaL)\n\n# Building CaN current\nDA_CaSCurrent = initializeCurrent(\"CaN\", DA_ECa, exponents=1,\n    activationSteadyStateGating=DA_mCaN_inf, activationTimeConstant=DA_tau_mCaN)\n\n# Building ERG current\nDA_ERGCurrent = initializeCurrent(\"ERG\", DA_EK, exponents=1,\n    activationSteadyStateGating=DA_o_inf, activationTimeConstant=DA_tau_o)\n\n# Building NMDA current\nDA_NMDACurrent = initializeCurrent(\"NMDA\", DA_ENMDA, exponents=1,\n    activationSteadyStateGating=DA_NMDA_inf, activationTimeConstant=DA_tau_NMDA,\n    MgDependency=true)\n\n# Building a more complex model with calcium\nDA_ionCurrents = [DA_NaCurrent, DA_KdCurrent, DA_CaLCurrent, DA_CaSCurrent, DA_ERGCurrent, DA_NMDACurrent]\nDA = initializeNeuronModel(DA_ionCurrents)","category":"page"},{"location":"initializing/examples/","page":"Examples of existing models","title":"Examples of existing models","text":"Note that here, we do not specify any membrane capacitance, nor leakage parameters or maximum ion channel conductances. Doing so will put default values for C, g_mathrmleak and E_mathrmleak, being respectively 1, 1 and -50. Not initializing any maximum ion channel conductances just indicate that the field maximumConductances of data structure NeuronCB will be filled with NaN. This means that no maximum ion channel conductances is specified for such model.","category":"page"},{"location":"#[NmodController.jl:-A-simple-way-to-use-conductance-based-models-and-to-simulate-them](https://github.com/arthur-fyon/NmodController.jl)","page":"Home","title":"NmodController.jl: A simple way to use conductance based models and to simulate them","text":"","category":"section"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"#Users","page":"Home","title":"Users","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Download Julia v1.6.X or later, if you haven't already.\nAdd the NmodController module entering the following at the Julia REPL ]add https://github.com/arthur-fyon/NmodController.jl.","category":"page"},{"location":"#Developers","page":"Home","title":"Developers","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Clone the NmodController module to username/.julia/dev/.\nEnter the package manager in REPL by pressing ] then add the package by typing dev NmodController rather than add NmodController.","category":"page"},{"location":"#Table-of-Contents","page":"Home","title":"Table of Contents","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Pages = [\"initializing/current.md\",\n         \"initializing/neuron.md\",\n         \"initializing/examples.md\",\n         \"computing/DICs.md\",\n         \"computing/examples.md\",\n         \"simulating/files.md\",\n         \"simulating/examples.md\",          \n         \"typesMethodsFunctions.md\"]","category":"page"},{"location":"initializing/neuron/#Initializing-a-conductance-based-model","page":"How to initialize a conductance based model?","title":"Initializing a conductance based model","text":"","category":"section"},{"location":"initializing/neuron/#Classical-conductance-based-model-without-intracellular-calcium-dynamic","page":"How to initialize a conductance based model?","title":"Classical conductance based model without intracellular calcium dynamic","text":"","category":"section"},{"location":"initializing/neuron/","page":"How to initialize a conductance based model?","title":"How to initialize a conductance based model?","text":"Once all the ionic current data structures have been initialized, any complete conductance based model can be descrbied. The voltage equation for such model, without any external applied current, writes","category":"page"},{"location":"initializing/neuron/","page":"How to initialize a conductance based model?","title":"How to initialize a conductance based model?","text":"dotV = (1C) cdot (-sum_mathrmion I_mathrmion - I_mathrmleak)","category":"page"},{"location":"initializing/neuron/","page":"How to initialize a conductance based model?","title":"How to initialize a conductance based model?","text":"where C is the membrane capacitance and I_mathrmleak writes","category":"page"},{"location":"initializing/neuron/","page":"How to initialize a conductance based model?","title":"How to initialize a conductance based model?","text":"I_mathrmleak = g_mathrmleak (V - E_mathrmleak)","category":"page"},{"location":"initializing/neuron/","page":"How to initialize a conductance based model?","title":"How to initialize a conductance based model?","text":"All the other equations in the model consists in the gating variables dynamics previously described.","category":"page"},{"location":"initializing/neuron/","page":"How to initialize a conductance based model?","title":"How to initialize a conductance based model?","text":"In NmodController.jl, a conductance based model can be contained in a NeuronCB type. To help initializing such data structure, calling initializeNeuronModel() with appropriate arguments is strongly recommended.","category":"page"},{"location":"initializing/neuron/#Example-1","page":"How to initialize a conductance based model?","title":"Example 1","text":"","category":"section"},{"location":"initializing/neuron/","page":"How to initialize a conductance based model?","title":"How to initialize a conductance based model?","text":"The next few lines of code show how to initialize a conductance based model with two ionic currents: a fast sodium and a rectified delayed potassium, i.e. the original Hodgkin and Huxley model. This voltage equation writes","category":"page"},{"location":"initializing/neuron/","page":"How to initialize a conductance based model?","title":"How to initialize a conductance based model?","text":"dotV = (1C) cdot (-I_mathrmNa -  I_mathrmKd - I_mathrmleak)","category":"page"},{"location":"initializing/neuron/","page":"How to initialize a conductance based model?","title":"How to initialize a conductance based model?","text":"# First wrapping all ionic currents in a vector\nionCurrents = [NaCurrent, KdCurrent]\n\n# Initializing leakage reversal potential and leakage conductance\nEleak = -50.\ngleak = 0.01\n\n# Initializing all the maximum ion channel conductances\nbar_g = [100., 10.]\n\nusing NmodController\nHHmodel = initializeNeuronModel(ionCurrents, C=0.1, leakageConductance=gleak, reversaleLeakagePotential=Eleak, maximumConductances=bar_g)","category":"page"},{"location":"initializing/neuron/","page":"How to initialize a conductance based model?","title":"How to initialize a conductance based model?","text":"Note that the argument maximumConductances is optional and correspond to all the maximum ion channel conductances barg_mathrmion wrapped in a vector. maximumConductances will be filled with NaN if not provided in input.","category":"page"},{"location":"initializing/neuron/#Intracellular-calcium-dynamic-in-a-conductance-based-model","page":"How to initialize a conductance based model?","title":"Intracellular calcium dynamic in a conductance based model","text":"","category":"section"},{"location":"initializing/neuron/","page":"How to initialize a conductance based model?","title":"How to initialize a conductance based model?","text":"When at least one ionic current is calcium dependent, an additional ODE has to be added to the conductance based model to describe intracellular calcium dynamics. Such equation generally writes","category":"page"},{"location":"initializing/neuron/","page":"How to initialize a conductance based model?","title":"How to initialize a conductance based model?","text":"tau_mathrmCa cdot dotCa = sum_mathrmionCa e_mathrmionCa I_mathrmionCa - Ca + Ca_mathrmeq","category":"page"},{"location":"initializing/neuron/","page":"How to initialize a conductance based model?","title":"How to initialize a conductance based model?","text":"where I_mathrmionCa is a calcium current of the model and e_mathrmionCa is its associated coefficient.","category":"page"},{"location":"initializing/neuron/","page":"How to initialize a conductance based model?","title":"How to initialize a conductance based model?","text":"In NmodController.jl, such intracellular calcium dynamics can be contained in a CalciumDynamic type. To help initializing such data structure, calling initializeCalciumDynamics() with appropriate arguments is strongly recommended.","category":"page"},{"location":"initializing/neuron/#Example-2","page":"How to initialize a conductance based model?","title":"Example 2","text":"","category":"section"},{"location":"initializing/neuron/","page":"How to initialize a conductance based model?","title":"How to initialize a conductance based model?","text":"The next few lines of code show how to initialize a conductance based model with five ionic currents: a fast sodium, a rectified delayed potassium, a calcium controlled potassium current, a T-type calcium current and a slow calcium current. Note that the calcium controlled potassium current is a special current that depends on both the voltage and intracellular calcium. This voltage equation writes","category":"page"},{"location":"initializing/neuron/","page":"How to initialize a conductance based model?","title":"How to initialize a conductance based model?","text":"dotV = (1C) cdot (-I_mathrmNa -  I_mathrmKd -  I_mathrmKCa -  I_mathrmCaS -  I_mathrmCaT - I_mathrmleak)","category":"page"},{"location":"initializing/neuron/","page":"How to initialize a conductance based model?","title":"How to initialize a conductance based model?","text":"The associated calcium dynamic writes","category":"page"},{"location":"initializing/neuron/","page":"How to initialize a conductance based model?","title":"How to initialize a conductance based model?","text":"tau_mathrmCa cdot dotCa = -094cdot  I_mathrmCaT -094cdot  I_mathrmCaS - Ca + Ca_mathrmeq","category":"page"},{"location":"initializing/neuron/","page":"How to initialize a conductance based model?","title":"How to initialize a conductance based model?","text":"with tau_mathrmCa = 20 and Ca_mathrmeq = 005.","category":"page"},{"location":"initializing/neuron/","page":"How to initialize a conductance based model?","title":"How to initialize a conductance based model?","text":"# First wrapping all ionic currents in a vector\nionCurrents = [NaCurrent, KdCurrent, KCaCurrent, CaTCurrent, CaSCurrent]\n\n# Initializing leakage reversal potential and leakage conductance\nEleak = -50.\ngleak = 0.01\n\n# Initializing all the maximum ion channel conductances\nbar_g = [100., 10., 10., 1., 1.]\n\nusing NmodController\nCaDyn = initializeCalciumDynamics([\"CaT\", \"CaS\"], [-0.94, -0.94], 0.05, 20)\nCBmodel = initializeNeuronModel(ionCurrents, C=0.1, calciumDynamics=CaDyn, leakageConductance=gleak, reversaleLeakagePotential=Eleak, maximumConductances=bar_g)","category":"page"}]
}
