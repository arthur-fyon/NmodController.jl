module NmodController
# Calling dependencies
using SymPy

# Calling all the files of the package
include("mathFunctions.jl")
include("neuron.jl")
include("elementsOfS.jl")
include("DIC.jl")
include("writeSimuFiles.jl")

# Export meaningful functions
export boltz
export tauX
export IonCurrent
export CalciumDynamic
export NeuronCB
export initializeCurrent
export initializeCalciumDynamics
export initializeNeuronModel
export computeDICs
export computeThresholdVoltage
export writeUncontrolledODEs
export writeControlledODEs

end
