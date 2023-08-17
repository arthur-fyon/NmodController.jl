module NmodController
# Calling dependencies
using SymPy

# Calling all the files of the package
include("neuron.jl")

# Export meaningful functions
export boltz
export tauX
export IonCurrent
export CalciumDynamic
export NeuronCB
export initializeCurrent
export initializeCalciumDynamics
export initializeNeuronModel

end
