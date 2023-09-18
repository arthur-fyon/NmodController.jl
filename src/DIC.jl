"""
    computeThresholdVoltage(gf::Function, 
        gs::Function, 
        gu::Function)

Compute threshold voltage of a specific conductance based model.
Refer to Drion et al., 2015 "Dynamic Input Conductances Shape Neuronal Spiking"

# Arguments
- `gf`: fast DIC. Must be of type Function.
- `gs`: slow DIC. Must be of type Function.
- `gu`: ultraslow DIC. Must be of type Function.

# Example
```jldoctest
julia> Vth = computeThresholdVoltage(gf, gs, gu)
-48.897088456154556
```
"""
function computeThresholdVoltage(gf::Function, gs::Function, gu::Function)
    # First compute the input conductance
    gin(V) = gf(V) + gs(V) + gu(V)

    # Initialize the two boundaries for the bisection algorithm
    lowBound = -80.
    highBound = lowBound

    # First check that gin is positive. If not, bring it to positive value
    if gin(lowBound) < 0.
        while gin(highBound) < 0.
            highBound += 0.1
    
            # If no positive value of gin, throw error
            if highBound > 80.
                error("No threshold voltage found for this model, no positive to negative crossing of input conductance!")
            end
        end

        # Put the low bound of the bisection where gin gets positive
        lowBound = highBound
    end

    # While the input conductance is still positive, go further
    while gin(highBound) > 0.
        highBound += 0.1

        # If no negative value of gin, throw error
        if highBound > 80.
            error("No threshold voltage found for this model, no positive to negative crossing of input conductance!")
        end
    end

    Vth = bisection(gin, lowBound, highBound)
    return Vth
end

"""
    computeDICs(neuron::NeuronCB, 
        tauFast::Function, 
        tauSlow::Function, 
        tauUltraslow::Function; 
        tauCa::Union{Float64, Int64, Function}=NaN, 
        Mg::Union{Int64, Float64}=NaN, 
        onlyS::Bool=false, 
        scaled::Bool=true)

Compute dynamic input conductances gDIC (DICs) or only the sensitivity matrix S of a conductance based model. Note that gDIC = S*bar_g with bar_g being the vector of all maximum ion channel conductances.
Refer to Drion et al., 2015 "Dynamic Input Conductances Shape Neuronal Spiking"

# Arguments
- `neuron`: data structure containing the conductance based model of interest.
- `tauFast`: fast time constant function, used to compute the DIC weights.
- `tauSlow`: slow time constant function, used to compute the DIC weights.
- `tauUltraslow`: ultraslow time constant function, used to compute the DIC weights.
- `tauCa`: time constant of the calcium if the model is calcium dependent, used to the DIC weights for the calcium only. Optional. Default value is NaN. Required if the neuron model is calcium dependent.
- `Mg`: magnesium concentration if the model is magnesium dependent. Optional. Default value is NaN. Required if the neuron model is magnesium dependent.
- `onlyS`: flag indicating if only the sensitivity matrix should be computed, for the case where no maximum ion channel conductances are in the model data structure. Optional. Default value is false.
- `scaled`: flag indicating if DICs or sensitivity matrix should be scaled by the leakage conductance. Optional. Default value is true.

# Example
```jldoctest
julia> gf, gs, gu = computeDICs(STG, STG_tauFast, STG_tauSlow, STG_tauUltraslow, tauCa=5000.);

julia> [gf(0), gs(0), gu(0)]
3-element Vector{Float64}:
  2712.186654147375
 19069.69205128652 
  -196.85745712674108

julia> S = computeDICs(STG, STG_tauFast, STG_tauSlow, STG_tauUltraslow, tauCa=5000., onlyS=true);

julia> S(0)
3×7 Matrix{Float64}:
 0.0213593  0.23756  0.00590276   0.000796009   0.0766643   29.8745   0.000857487
 0.059507   2.47411  0.0353302   -0.00824027    0.118022   211.202   -0.0
 0.0        1.29588  0.0386544   -0.00383275   -3.34257      0.0     -0.00285826
```
"""
function computeDICs(neuron::NeuronCB, tauFast::Function, tauSlow::Function, tauUltraslow::Function; tauCa::Union{Float64, Int64, Function}=NaN, Mg::Union{Int64, Float64}=NaN, 
    onlyS::Bool=false, scaled::Bool=true)
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

        if scaled
            # Scaling by leakage leakage conductance
            S_scaled(V) = S(V) ./ neuron.leakageConductance
            return S_scaled
        else
            return S
        end
    # Otherwise, compute the DICs
    else
        # Compute the submatrices of the sensitivity matrix
        Sf(V) = [S_Tuple[i][1](V) for i ∈ 1 : neuron.numberOfIonCurrents]
        Ss(V) = [S_Tuple[i][2](V) for i ∈ 1 : neuron.numberOfIonCurrents]
        Su(V) = [S_Tuple[i][3](V) for i ∈ 1 : neuron.numberOfIonCurrents]

        if scaled
            # Computing normalized DICs
            gf_scaled(V) = Sf(V)' * neuron.maximumConductances ./ neuron.leakageConductance  + 1.
            gs_scaled(V) = Ss(V)' * neuron.maximumConductances ./ neuron.leakageConductance
            gu_scaled(V) = Su(V)' * neuron.maximumConductances ./ neuron.leakageConductance

            return gf_scaled, gs_scaled, gu_scaled
        else
            # Computing unnormalized DICs
            gf(V) = Sf(V)' * neuron.maximumConductances  + neuron.leakageConductance
            gs(V) = Ss(V)' * neuron.maximumConductances
            gu(V) = Su(V)' * neuron.maximumConductances

            return gf, gs, gu
        end
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
