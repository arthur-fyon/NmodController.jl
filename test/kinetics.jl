### Kinetics and parameters for the STG model
# Initializing Nernst reversalPotential
STG_VNa = 50. # Sodium reversal potential
STG_VK = -80. # Potassium reversal potential
STG_VCa = 80. # Calcium reversal potential
STG_VH = -20. # Reversal potential for the H-current (permeable to both sodium and potassium ions)
STG_Vleak = -50. # Reversal potential of leak channels

# Na current
STG_mNa_inf(V) = boltz(V, 25.5, -5.29)
STG_tau_mNa(V) = tauX(V, 1.32, 1.26, 120., -25.)
STG_hNa_inf(V) = boltz(V, 48.9, 5.18)
STG_tau_hNa(V) = (0.67 / (1 + exp((V+62.9)/-10.0))) * (1.5 + 1 / (1+exp((V+34.9)/3.6)))

# Kd current
STG_mKd_inf(V) = boltz(V, 12.3, -11.8)
STG_tau_mKd(V) = tauX(V, 7.2, 6.4, 28.3, -19.2)

# KCa current
STG_mKCa_inf(V, Ca) = (Ca / (Ca+3)) * (1 / (1+exp((V+28.3)/-12.6)))
STG_tau_mKCa(V) = tauX(V, 90.3, 75.1, 46., -22.7)

# CaT current
STG_mCaT_inf(V) = boltz(V, 27.1, -7.2)
STG_tau_mCaT(V) = tauX(V, 21.7, 21.3, 68.1, -20.5)
STG_hCaT_inf(V) = boltz(V, 32.1, 5.5)
STG_tau_hCaT(V) = tauX(V, 105., 89.8, 55., -16.9)

# CaS current
STG_mCaS_inf(V) = boltz(V, 33., -8.1)
STG_tau_mCaS(V) = 1.4 + (7 / ((exp((V+27)/10)) + (exp((V+70)/-13))))
STG_hCaS_inf(V) = boltz(V, 60., 6.2)
STG_tau_hCaS(V) = 60 + (150 / ((exp((V+55)/9)) + (exp((V+65)/-16))))

# A current
STG_mA_inf(V) = boltz(V, 27.2, -8.7)
STG_tau_mA(V) = tauX(V, 11.6, 10.4, 32.9, -15.2)
STG_hA_inf(V) = boltz(V, 56.9, 4.9)
STG_tau_hA(V) = tauX(V, 38.6, 29.2, 38.9, -26.5)

# H current
STG_mH_inf(V) = boltz(V, 70., 6.)
STG_tau_mH(V) = tauX(V, 272., -1499., 42.2, -8.73)

### Kinetics and parameters for the DA model
# Initializing Nernst reversalPotential
DA_VNa = 60. # Sodium reversal potential
DA_VK = -85. # Potassium reversal potential
DA_VCa = 60. # Calcium reversal potential
DA_VNMDA = 0. # NMDA reversal potential
DA_Vleak = -50. # Reversal potential of leak channels
DA_Mg = 1.4 # Mg concentration

# Na current
DA_mNa_inf(V) = boltz(V, -30.0907, 9.7264)
DA_tau_mNa(V) = 0.01 + 1.0 / ((-(15.6504 + 0.4043*V)/(exp(-19.565 -0.5052*V)-1.0)) + 3.0212*exp(-7.4630e-3*V))
DA_hNa_inf(V) = boltz(V, -54.0289, -10.7665)
DA_tau_hNa(V) = 0.4 + 1.0 / ((5.0754e-4*exp(-6.3213e-2*V)) + 9.7529*exp(0.13442*V))

# Kd current
DA_mKd_inf(V) = boltz(V, -25., 12.)
DA_tau_mKd(V) = tauX(V, 20., 18., 38., -10.)

# CaL current
DA_mCaL_inf(V) = boltz(V, -50., 2.)
DA_tau_mCaL(V) = tauX(V, 30., 28., 45., -3.)

# CaN current
DA_mCaN_inf(V) = boltz(V, -30., 7.)
DA_tau_mCaN(V) = tauX(V, 30., 25., 55., -6.)

# ERG current
a0ERG(V) = 0.0036 * exp(0.0759*V)
b0ERG(V) = 1.2523e-5 * exp(-0.0671*V)
aiERG(V) = 0.1 * exp(0.1189*V)
biERG(V) = 0.003 * exp(-0.0733*V)
DA_o_inf(V) = a0ERG(V)*biERG(V) / (a0ERG(V)*(aiERG(V)+biERG(V)) + b0ERG(V)*biERG(V))
DA_tau_o(V) = 100000.

# NMDA current
DA_NMDA_inf(V, Mg) = 1 / (1 + Mg*exp(-0.08*V)/10.)
DA_tau_NMDA(V) = 1e-10
