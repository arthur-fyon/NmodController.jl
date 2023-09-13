var documenterSearchIndex = {"docs":
[{"location":"computing/DICs/#WIP","page":"How to compute DICs and sensitivity matrix?","title":"WIP","text":"","category":"section"},{"location":"initializing/current/#Initializing-an-ionic-current","page":"How to initialize an ionic current?","title":"Initializing an ionic current","text":"","category":"section"},{"location":"initializing/current/","page":"How to initialize an ionic current?","title":"How to initialize an ionic current?","text":"In conductance based models, an ionic current I_mathrmion is described by its voltage dependent conductance g_mathrmion and Nernst reversal potential E_mathrmion","category":"page"},{"location":"initializing/current/","page":"How to initialize an ionic current?","title":"How to initialize an ionic current?","text":"I_mathrmion = g_mathrmion(V) cdot (V - E_mathrmion)","category":"page"},{"location":"initializing/current/","page":"How to initialize an ionic current?","title":"How to initialize an ionic current?","text":"This voltage dependent conductance is defined by its maximum value barg_mathrmion and at most two gating variables, one activation m_mathrmion(V) and one inactivation h_mathrmion(V), that varies between 0 and 1. In the case of a current with two gating variables, this writes","category":"page"},{"location":"initializing/current/","page":"How to initialize an ionic current?","title":"How to initialize an ionic current?","text":"g_mathrmion(V) = barg_mathrmion cdot m^a_mathrmion(V) cdot h^b_mathrmion(V)","category":"page"},{"location":"initializing/current/","page":"How to initialize an ionic current?","title":"How to initialize an ionic current?","text":"Note that the maximum ion channel conductance of the current will not be stocked inside the current data structure, but rather in the conductance based model data structure.","category":"page"},{"location":"initializing/current/","page":"How to initialize an ionic current?","title":"How to initialize an ionic current?","text":"Each gating variable dynamic follow a basic first order ODE where both the time constant tau_m_mathrmion and the converging value m_mathrmion_infty are voltage dependent: ","category":"page"},{"location":"initializing/current/","page":"How to initialize an ionic current?","title":"How to initialize an ionic current?","text":"dotm_mathrmion = fracm_mathrmion_infty(V) - m_mathrmiontau_m_mathrmion(V)","category":"page"},{"location":"initializing/current/","page":"How to initialize an ionic current?","title":"How to initialize an ionic current?","text":"In NmodController.jl, an ionic current can be contained in a IonCurrent type. To help initializing such data structure, calling initializeCurrent() with appropriate arguments is strongly recommended.","category":"page"},{"location":"initializing/current/#Example-1","page":"How to initialize an ionic current?","title":"Example 1","text":"","category":"section"},{"location":"initializing/current/","page":"How to initialize an ionic current?","title":"How to initialize an ionic current?","text":"The next few lines of code show how to initialize a sodium current with two gating variables. This current writes","category":"page"},{"location":"initializing/current/","page":"How to initialize an ionic current?","title":"How to initialize an ionic current?","text":"I_mathrmNa = barg_mathrmNa cdot m^3_mathrmNa cdot h_mathrmNa cdot (V - E_mathrmNa)","category":"page"},{"location":"initializing/current/","page":"How to initialize an ionic current?","title":"How to initialize an ionic current?","text":"# First initializing the converging values and time constants functions\nboltz(V, A, B) = 1 / (1 + exp((V+A) / B))\ntauX(V, A, B, D, E) = A - B / (1 + exp((V+D) / E))\n\nmNa_inf(V) = boltz(V, 25.5, -5.29)\ntau_mNa(V) = tauX(V, 1.32, 1.26, 120., -25.)\nhNa_inf(V) = boltz(V, 48.9, 5.18)\ntau_hNa(V) = (0.67 / (1 + exp((V+62.9)/-10.0))) * (1.5 + 1 / (1+exp((V+34.9)/3.6)))\n\n# Initializing Nernst reversal potential\nENa = 50.\n\nusing NmodController\nNaCurrent = initializeCurrent(\"Na\", ENa, numberOfGatings=2, exponents=[3, 1],\n    activationSteadyStateGating=mNa_inf, activationTimeConstant=tau_mNa,\n    inactivationSteadyStateGating=hNa_inf, inactivationTimeConstant=tau_hNa)","category":"page"},{"location":"initializing/current/#Example-2","page":"How to initialize an ionic current?","title":"Example 2","text":"","category":"section"},{"location":"initializing/current/","page":"How to initialize an ionic current?","title":"How to initialize an ionic current?","text":"Sometimes, the gating variable converging values may depend on the intracellular calcium. The next few lines of code show how to initialize a calcium controlled potassium current with one gating variable. This current writes","category":"page"},{"location":"initializing/current/","page":"How to initialize an ionic current?","title":"How to initialize an ionic current?","text":"I_mathrmKCa = barg_mathrmKCa cdot m^4_mathrmKCa cdot (V - E_mathrmKCa)","category":"page"},{"location":"initializing/current/","page":"How to initialize an ionic current?","title":"How to initialize an ionic current?","text":"Where the dynamic of m_mathrmKCa is described by","category":"page"},{"location":"initializing/current/","page":"How to initialize an ionic current?","title":"How to initialize an ionic current?","text":"dotm_mathrmKCa = fracm_mathrmKCa_infty(V Ca) - m_mathrmKCatau_m_mathrmKCa(V)","category":"page"},{"location":"initializing/current/","page":"How to initialize an ionic current?","title":"How to initialize an ionic current?","text":"As you can see, m_mathrmKCa_infty depends on both the voltage and the calcium, this must be notified to initializeCurrent().","category":"page"},{"location":"initializing/current/","page":"How to initialize an ionic current?","title":"How to initialize an ionic current?","text":"# First initializing the converging value and time constant functions\nboltz(V, A, B) = 1 / (1 + exp((V+A) / B))\ntauX(V, A, B, D, E) = A - B / (1 + exp((V+D) / E))\n\nmKCa_inf(V, Ca) = (Ca / (Ca+3)) * (1 / (1+exp((V+28.3)/-12.6)))\ntau_mKCa(V) = tauX(V, 90.3, 75.1, 46., -22.7)\n\n# Initializing Nernst reversal potential\nEK = -80.\n\nusing NmodController\nKCaCurrent = initializeCurrent(\"KCa\", EK, exponents=4,\n    activationSteadyStateGating=mKCa_inf, activationTimeConstant=tau_mKCa,\n    calciumDependency=true)","category":"page"},{"location":"initializing/current/#Example-3","page":"How to initialize an ionic current?","title":"Example 3","text":"","category":"section"},{"location":"initializing/current/","page":"How to initialize an ionic current?","title":"How to initialize an ionic current?","text":"Another specific case is where the time constant may not depend on the voltage. In such case, just provide the time constant as a Float64 or a Int64 in argument of initializeCurrent(). Moreover, the converging value of the gating variable may be magnesium dependent. In such case, do as the following. The next few lines of code show how to initialize an instantaneous magnesium dependent NMDA current with one gating variable. This current writes","category":"page"},{"location":"initializing/current/","page":"How to initialize an ionic current?","title":"How to initialize an ionic current?","text":"I_mathrmNMDA = barg_mathrmNMDA cdot m_mathrmNMDA cdot (V - E_mathrmNMDA)","category":"page"},{"location":"initializing/current/","page":"How to initialize an ionic current?","title":"How to initialize an ionic current?","text":"Where the dynamic of m_mathrmNMDA is described by","category":"page"},{"location":"initializing/current/","page":"How to initialize an ionic current?","title":"How to initialize an ionic current?","text":"m_mathrmNMDA = m_mathrmNMDA_infty(V Mg)","category":"page"},{"location":"initializing/current/","page":"How to initialize an ionic current?","title":"How to initialize an ionic current?","text":"As you can see, m_mathrmNMDA_infty depends on both the voltage and the magnesium, this must be notified to initializeCurrent()","category":"page"},{"location":"initializing/current/","page":"How to initialize an ionic current?","title":"How to initialize an ionic current?","text":"# First initializing the converging value and time constant functions\nNMDA_inf(V, Mg) = 1 / (1 + Mg*exp(-0.08*V)/10.)\ntau_NMDA(V) = 1e-10 # Never equal to zero, this might be fixed later on\n\n# Initializing Nernst reversal potential\nENMDA = 0.\n\nusing NmodController\nNMDACurrent = initializeCurrent(\"NMDA\", ENMDA, exponents=1,\n    activationSteadyStateGating=NMDA_inf, activationTimeConstant=tau_NMDA,\n    MgDependency=true)","category":"page"},{"location":"computing/examples/#WIP","page":"Examples of existing models","title":"WIP","text":"","category":"section"},{"location":"simulating/examples/#WIP","page":"Examples of existing models","title":"WIP","text":"","category":"section"},{"location":"simulating/files/#WIP","page":"How to write ODE functions file?","title":"WIP","text":"","category":"section"},{"location":"initializing/examples/#Example-of-existing-models-in-the-litterature","page":"Examples of existing models","title":"Example of existing models in the litterature","text":"","category":"section"},{"location":"initializing/examples/","page":"Examples of existing models","title":"Examples of existing models","text":"Through this documentation, two conductance based model coming from the litterature wil be used with NmodController.jl: a stomatogastric (STG) neuron model (Liu et al., 1998 \"A model neuron with activity-dependent conductances regulated by multiple calcium sensors\") and an adapted dopaminergic (DA) neuron model (Drion et al., 2011 \"How modeling can reconcile apparently discrepant experimental results: the case of pacemaking in dopaminergic neurons\").","category":"page"},{"location":"initializing/examples/#Initializing-the-STG-model","page":"Examples of existing models","title":"Initializing the STG model","text":"","category":"section"},{"location":"initializing/examples/","page":"Examples of existing models","title":"Examples of existing models","text":"The STG model is composed of 7 voltage gated ionic currents which one is calcium dependent:","category":"page"},{"location":"initializing/examples/","page":"Examples of existing models","title":"Examples of existing models","text":"Transient sodium current I_mathrmNa (2 gating variables);\nT-type calcium current I_mathrmCaT (2 gating variables);\nSlow calcium current I_mathrmCaS (2 gating variables);\nA-type potassium current I_mathrmA (2 gating variables);\nCalcium controlled potassium current I_mathrmKCa (1 gating variable caclium dpendent);\nDelayed rectified potassium current I_mathrmKd (1 gating variable);\nH current I_mathrmH (1 gating variable).","category":"page"},{"location":"initializing/examples/","page":"Examples of existing models","title":"Examples of existing models","text":"The voltage equation writes","category":"page"},{"location":"initializing/examples/","page":"Examples of existing models","title":"Examples of existing models","text":"$","category":"page"},{"location":"initializing/examples/","page":"Examples of existing models","title":"Examples of existing models","text":"\\begin{eqnarray} C \\dot V = &-& \\bar{g}\\mathrm{Na}m^3\\mathrm{Na}h\\mathrm{Na}(V-E\\mathrm{Na}) \\             &-& \\bar{g}\\mathrm{CaT}m^3\\mathrm{CaT}h\\mathrm{CaT}(V-E\\mathrm{Ca}) \\             &-& \\bar{g}\\mathrm{CaS}m^3\\mathrm{CaS}h\\mathrm{CaS}(V-E\\mathrm{Ca}) \\             &-& \\bar{g}\\mathrm{A}m^3\\mathrm{A}h\\mathrm{A}(V-E\\mathrm{K}) \\             &-& \\bar{g}\\mathrm{KCa}m^4\\mathrm{KCa}(V-E\\mathrm{K}) \\             &-& \\bar{g}\\mathrm{Kd}m^4\\mathrm{Kd}(V-E\\mathrm{K}) \\             &-& \\bar{g}\\mathrm{H}m\\mathrm{H}(V-E\\mathrm{H}) \\             &-& g\\mathrm{leak}(V-E\\mathrm{leak}) \\\n           &+& I{ext}(t) \\end{eqnarray} $","category":"page"},{"location":"initializing/examples/","page":"Examples of existing models","title":"Examples of existing models","text":"and the intracellular calcium dynamic writes","category":"page"},{"location":"initializing/examples/","page":"Examples of existing models","title":"Examples of existing models","text":"tau_Ca cdot dotCa = -094cdot  I_CaT -094cdot  I_CaS - Ca + Ca_mathrmeq","category":"page"},{"location":"initializing/examples/","page":"Examples of existing models","title":"Examples of existing models","text":"with tau_Ca = 20 and Ca_mathrmeq = 005.","category":"page"},{"location":"initializing/examples/","page":"Examples of existing models","title":"Examples of existing models","text":"The next few lines of code show how to initialize such model using NmodController.jl.","category":"page"},{"location":"#[NmodController.jl:-A-simple-way-to-use-conductance-based-models-and-to-simulate-them](https://github.com/arthur-fyon/NmodController.jl)","page":"Home","title":"NmodController.jl: A simple way to use conductance based models and to simulate them","text":"","category":"section"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"#Users","page":"Home","title":"Users","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Download Julia v1.6.X or later, if you haven't already.\nAdd the NmodController module entering the following at the Julia REPL ]add https://github.com/arthur-fyon/NmodController.jl.","category":"page"},{"location":"#Developers","page":"Home","title":"Developers","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Clone the NmodController module to username/.julia/dev/.\nEnter the package manager in REPL by pressing ] then add the package by typing dev NmodController rather than add NmodController.","category":"page"},{"location":"#Table-of-Contents","page":"Home","title":"Table of Contents","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Pages = [\"initializing/current.md\",\n         \"initializing/neuron.md\",\n         \"initializing/examples.md\",\n         \"computing/DICs.md\",\n         \"computing/examples.md\",\n         \"simulating/files.md\",\n         \"simulating/examples.md\",          \n         \"typesMethodsFunctions.md\"]","category":"page"},{"location":"initializing/neuron/#Initializing-a-conductance-based-model","page":"How to initialize a conductance based model?","title":"Initializing a conductance based model","text":"","category":"section"},{"location":"initializing/neuron/","page":"How to initialize a conductance based model?","title":"How to initialize a conductance based model?","text":"Once all the ionic current data structures have been initialized, any complete conductance based model can be descrbied. The voltage equation for such model, without any external applied current, writes","category":"page"},{"location":"initializing/neuron/","page":"How to initialize a conductance based model?","title":"How to initialize a conductance based model?","text":"dotV = (1C) cdot (-sum_mathrmion I_mathrmion - I_mathrmleak)","category":"page"},{"location":"initializing/neuron/","page":"How to initialize a conductance based model?","title":"How to initialize a conductance based model?","text":"where C is the membrane capacitance and I_mathrmleak writes","category":"page"},{"location":"initializing/neuron/","page":"How to initialize a conductance based model?","title":"How to initialize a conductance based model?","text":"I_mathrmleak = g_mathrmleak (V - E_mathrmleak)","category":"page"},{"location":"initializing/neuron/","page":"How to initialize a conductance based model?","title":"How to initialize a conductance based model?","text":"All the other equations in the model consists in the gating variables dynamics previously described.","category":"page"},{"location":"initializing/neuron/","page":"How to initialize a conductance based model?","title":"How to initialize a conductance based model?","text":"In NmodController.jl, a conductance based model can be contained in a NeuronCB type. To help initializing such data structure, calling initializeNeuronModel() with appropriate arguments is strongly recommended.","category":"page"},{"location":"initializing/neuron/#Example-1","page":"How to initialize a conductance based model?","title":"Example 1","text":"","category":"section"},{"location":"initializing/neuron/","page":"How to initialize a conductance based model?","title":"How to initialize a conductance based model?","text":"The next few lines of code show how to initialize a conductance based model with two ionic currents: a fast sodium and a rectified delayed potassium, i.e. the original Hodgkin and Huxley model. This voltage equation writes","category":"page"},{"location":"initializing/neuron/","page":"How to initialize a conductance based model?","title":"How to initialize a conductance based model?","text":"dotV = (1C) cdot (-I_mathrmNa -  I_mathrmKd - I_mathrmleak)","category":"page"},{"location":"initializing/neuron/","page":"How to initialize a conductance based model?","title":"How to initialize a conductance based model?","text":"# First wrapping all ionic currents in a vector\nionCurrents = [NaCurrent, KdCurrent]\n\n# Initializing leakage reversal potential and leakage conductance\nEleak = -50.\ngleak = 0.01\n\n# Initializing all the maximum ion channel conductances\nbar_g = [100., 10.]\n\nusing NmodController\nHHmodel = initializeNeuronModel(ionCurrents, C=0.1, leakageConductance=gleak, reversaleLeakagePotential=Eleak, maximumConductances=bar_g)","category":"page"},{"location":"initializing/neuron/","page":"How to initialize a conductance based model?","title":"How to initialize a conductance based model?","text":"Note that the argument maximumConductances is optional and correspond to all the maximum ion channel conductances barg_mathrmion wrapped in a vector. maximumConductances will be filled with NaN if not provided in input.","category":"page"},{"location":"initializing/neuron/#Intracellular-calcium-dynamic-in-a-conductance-based-model","page":"How to initialize a conductance based model?","title":"Intracellular calcium dynamic in a conductance based model","text":"","category":"section"},{"location":"initializing/neuron/","page":"How to initialize a conductance based model?","title":"How to initialize a conductance based model?","text":"When at least one ionic current is calcium dependent, an additional ODE has to be added to the conductance based model to describe intracellular calcium dynamics. Such equation generally writes","category":"page"},{"location":"initializing/neuron/","page":"How to initialize a conductance based model?","title":"How to initialize a conductance based model?","text":"tau_mathrmCa cdot dotCa = sum_mathrmionCa e_mathrmionCa I_mathrmionCa - Ca + Ca_mathrmeq","category":"page"},{"location":"initializing/neuron/","page":"How to initialize a conductance based model?","title":"How to initialize a conductance based model?","text":"where I_mathrmionCa is a calcium current of the model and e_mathrmionCa is its associated coefficient.","category":"page"},{"location":"initializing/neuron/","page":"How to initialize a conductance based model?","title":"How to initialize a conductance based model?","text":"In NmodController.jl, such intracellular calcium dynamics can be contained in a CalciumDynamic type. To help initializing such data structure, calling initializeCalciumDynamics() with appropriate arguments is strongly recommended.","category":"page"},{"location":"initializing/neuron/#Example-2","page":"How to initialize a conductance based model?","title":"Example 2","text":"","category":"section"},{"location":"initializing/neuron/","page":"How to initialize a conductance based model?","title":"How to initialize a conductance based model?","text":"The next few lines of code show how to initialize a conductance based model with five ionic currents: a fast sodium, a rectified delayed potassium, a calcium controlled potassium current, a T-type calcium current and a slow calcium current. Note that the calcium controlled potassium current is a special current that depends on both the voltage and intracellular calcium. This voltage equation writes","category":"page"},{"location":"initializing/neuron/","page":"How to initialize a conductance based model?","title":"How to initialize a conductance based model?","text":"dotV = (1C) cdot (-I_mathrmNa -  I_mathrmKd -  I_mathrmKCa -  I_mathrmCaS -  I_mathrmCaT - I_mathrmleak)","category":"page"},{"location":"initializing/neuron/","page":"How to initialize a conductance based model?","title":"How to initialize a conductance based model?","text":"The associated calcium dynamic writes","category":"page"},{"location":"initializing/neuron/","page":"How to initialize a conductance based model?","title":"How to initialize a conductance based model?","text":"tau_mathrmCa cdot dotCa = -094cdot  I_mathrmCaT -094cdot  I_mathrmCaS - Ca + Ca_mathrmeq","category":"page"},{"location":"initializing/neuron/","page":"How to initialize a conductance based model?","title":"How to initialize a conductance based model?","text":"with tau_mathrmCa = 20 and Ca_mathrmeq = 005.","category":"page"},{"location":"initializing/neuron/","page":"How to initialize a conductance based model?","title":"How to initialize a conductance based model?","text":"# First wrapping all ionic currents in a vector\nionCurrents = [NaCurrent, KdCurrent, KCaCurrent, CaTCurrent, CaSCurrent]\n\n# Initializing leakage reversal potential and leakage conductance\nEleak = -50.\ngleak = 0.01\n\n# Initializing all the maximum ion channel conductances\nbar_g = [100., 10., 10., 1., 1.]\n\n# Initialize the current dynamic data structure\nCaDyn = initializeCalciumDynamics([\"CaT\", \"CaS\"], [-0.94, -0.94], 0.05, 20)\n\nusing NmodController\nCBmodel = initializeNeuronModel(ionCurrents, C=0.1, calciumDynamics=CaDyn, leakageConductance=gleak, reversaleLeakagePotential=Eleak, maximumConductances=bar_g)","category":"page"}]
}
