# Function that computes the threshold voltage of the neuron
function computeThresholdVoltage(gf::Function, gs::Function, gu::Function)
    # First compute the input conductance
    gin(V) = gf(V) + gs(V) + gu(V)

    # Initialize the two boundaries for the bisection algorithm
    lowBound = -80.
    highBound = lowBound

    # While the input conductance is still positive, go further
    while gin(highBound) > 0.
        highBound += 0.1
    end

    Vth = bisection(gin, lowBound, highBound)
    return Vth
end

# Function that computes either sensitivity matrix or DICs
function computeDICs(neuron::NeuronCB, tauFast::Function, tauSlow::Function, tauUltraslow::Function; tauCa::Union{Float64, Int64, Function}=NaN, Mg::Union{Int64, Float64}=NaN, onlyS::Bool=false)
    # The neuron has not defined maximal conductances if the DICs have to be computed (onlyS==false)
    if !onlyS && any(isnan, neuron.maximumConductances)
        error("Please input a neuron which has defined values of maximal conductances, or compute only the sensitivity matrix by putting onlyS to true!")
    end

    # First compute the sensitivity matrix
    S_Tuple = computeS(neuron, tauFast, tauSlow, tauUltraslow, tauCa, Mg)

    # If only the sensitivity matrix is computed
    if onlyS
        # Reorganize the sensitivity matrix and make it callable
        S(V) = [S_Tuple[j][i](V) for i ∈ 1 : 3, j ∈ 1 : neuron.numberOfIonCurrents]

        # Scaling by leakage leakage conductance
        S_scaled(V) = S(V) ./ neuron.leakageConductance
        return S_scaled
    # Otherwise, compute the DICs
    else
        # Compute the submatrices of the sensitivity matrix and normalize them
        Sf(V) = [S_Tuple[i][1](V) for i ∈ 1 : neuron.numberOfIonCurrents] ./ neuron.leakageConductance
        Ss(V) = [S_Tuple[i][2](V) for i ∈ 1 : neuron.numberOfIonCurrents] ./ neuron.leakageConductance
        Su(V) = [S_Tuple[i][3](V) for i ∈ 1 : neuron.numberOfIonCurrents] ./ neuron.leakageConductance

        # Computing DICs
        gf(V) = Sf(V)' * neuron.maximumConductances + 1.
        gs(V) = Ss(V)' * neuron.maximumConductances
        gu(V) = Su(V)' * neuron.maximumConductances
        return gf, gs, gu
    end
end

# Function that computes the non normalized sensitivity matrix
function computeS(neuron::NeuronCB, tauFast::Function, tauSlow::Function, tauUltraslow::Function, tauCa::Union{Float64, Int64, Function}, Mg::Union{Int64, Float64})
    # First compute equilibrium value of calcium (0 if no calcium dependency in the model)
    CaEquilibrium = computeEquilibriumOfCalcium(neuron)
    
    # Compute the sensitivity matrix as a Vector of Tuple of Functions (fast, slow and ultraslow)
    S_Tuple = [computeWeightedElementOfS(ionCurrent, tauFast, tauSlow, tauUltraslow, CaEquilibrium, tauCa, Float64(Mg)) for ionCurrent in neuron.ionCurrents]

    return S_Tuple
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
