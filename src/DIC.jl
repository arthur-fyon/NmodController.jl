# Function that computes the DICs or the sensitivity matrix
function computeS(neuron::NeuronCB, tauFast::Function, tauSlow::Function, tauUltraslow::Function; tauCa::Union{Float64, Int64, Function}=NaN, Mg::Union{Int64, Float64}=NaN)
    # First compute equilibrium value of calcium (0 if no calcium dependency in the model)
    CaEquilibrium = computeEquilibriumOfCalcium(neuron)
    
    # Compute the sensitivity matrix as a Vector of Tuple of Functions (fast, slow and ultraslow)
    S_Tuple = [computeWeightedElementOfS(ionCurrent, tauFast, tauSlow, tauUltraslow, CaEquilibrium, tauCa, Float64(Mg)) for ionCurrent in neuron.ionCurrents]

    # Reorganize the sensitivity matrix and make it callable
    S(V) = [S_Tuple[j][i](V) for i ∈ 1 : 3, j ∈ 1 : neuron.numberOfIonCurrents]

    # Scaling by leakage leakage conductance
    S_scaled(V) = S(V) ./ neuron.leakageConductance
    return S_scaled
end

# Function that computes the DIC weights of one gating variable
function computeWeights(tauX::Function, tauFast::Function, tauSlow::Function, tauUltraslow::Function)
    # Compute logarithmic distances between adjacent time scales
    distFastSlow = logDist(tauX, tauFast, tauSlow)
    distSlowUslow  = logDist(tauX, tauSlow, tauUltraslow)
    
    # wfs = 1 if (tauX < tauFast) OR wfs = log distance if (tauFast <= tauX < tauSlow), 0 otherwise
    wfs(V) = 1. * (tauX(V) < tauFast(V)) + distFastSlow(V) * (tauFast(V) <= tauX(V) < tauSlow(V))

    # wsu = 1 if (tauX < tauSlow) OR wfs = log distance if (tauSlow <= tauX < tauUltraslow), 0 otherwise
    wsu(V) = 1. * (tauX(V) < tauSlow(V)) + distSlowUslow(V) * (tauSlow(V) <= tauX(V) < tauUltraslow(V))

    return wfs, wsu
end

# Function that computes the logarithmic distance
function logDist(tauX::Function, tauA::Function, tauB::Function)
    dist(V) = (log(tauB(V)) - log(tauX(V)))/(log(tauB(V)) - log(tauA(V)))
    return dist
end
