# Example on existing models in the litterature

# Computing DICs and threshold voltage on the STG model
Once the STG model is initialized, it is very straightforward to compute DICs, sensitivity matrix or threshold voltage. For comparison, two STG models with different maximum ion channel conductances will be used. `STG_spiking` that exhibits a tonic spiking behavior with
- sodium current: $\bar{g}_\mathrm{Na} =$ 4000;
- T-type calcium current: $\bar{g}_\mathrm{CaT} =$ 3;
- slow calcium current: $\bar{g}_\mathrm{CaS} =$ 4;
- A-type potassium current: $\bar{g}_\mathrm{A} =$ 175;
- calcium controlled potassium current: $\bar{g}_\mathrm{KCa} =$ 110;
- delayed rectified potassium current: $\bar{g}_\mathrm{Kd} =$ 137;
- H type current: $\bar{g}_\mathrm{H} =$ 0.3;
- leakage current: $g_\mathrm{leak} =$ 0.01;

and `STG_bursting` that exhibits a bursting behavior with
- sodium current: $\bar{g}_\mathrm{Na} =$ 4000;
- T-type calcium current: $\bar{g}_\mathrm{CaT} =$ 3;
- slow calcium current: $\bar{g}_\mathrm{CaS} =$ 19;
- A-type potassium current: $\bar{g}_\mathrm{A} =$ 70;
- calcium controlled potassium current: $\bar{g}_\mathrm{KCa} =$ 110;
- delayed rectified potassium current: $\bar{g}_\mathrm{Kd} =$ 137;
- H type current: $\bar{g}_\mathrm{H} =$ 0.3;
- leakage current: $g_\mathrm{leak} =$ 0.01.

The next few lines of code show how to compute and plot DICs, sensitivity matrix and threshold voltage of the spiking STG model using *NmodController.jl*.

```julia
# Defining some timescales
STG_tauFast = STG_tau_mNa
STG_tauSlow = STG_tau_mKd
STG_tauUltraslow = STG_tau_mH

# Computing DICs for the spiking model
STG_gf_spiking, STG_gs_spiking, STG_gu_spiking = computeDICs(STG_spiking, STG_tauFast, 
    STG_tauSlow, STG_tauUltraslow, tauCa=500000.)

# Computing the sensitivity matrix (which is a voltage function that can be called by STG_S_spiking(V))
STG_S_spiking = computeDICs(STG_spiking, STG_tauFast, 
    STG_tauSlow, STG_tauUltraslow, tauCa=500000., onlyS=true)

# Computing the threshold voltage
STG_Vth = computeThresholdVoltage(STG_gf_spiking, STG_gs_spiking, STG_gu_spiking)

# Plotting
V = -80 : 0.1 : 40.
Vzoom = -53. : 0.01 : -48.
p1 = plot(V, STG_gf_spiking.(V), linewidth=1.5, legend=false, size=(600, 200), color=:gray30)
vline!([Vzoom[1]], linewidth=1.5, color=:black, linestyle=:dash)
vline!([Vzoom[end]], linewidth=1.5, color=:black, linestyle=:dash)
ylabel!("g_f")
p1zoom = plot(Vzoom, STG_gf_spiking.(Vzoom), linewidth=1.5, legend=false, size=(300, 200), color=:gray30)
vline!([STG_Vth], linewidth=1.5, color=:firebrick1, linestyle=:dashdot)
p2 = plot(V, STG_gs_spiking.(V), linewidth=1.5, legend=false, size=(600, 200), color=:gray30)
vline!([Vzoom[1]], linewidth=1.5, color=:black, linestyle=:dash)
vline!([Vzoom[end]], linewidth=1.5, color=:black, linestyle=:dash)
ylabel!("g_s")
p2zoom = plot(Vzoom, STG_gs_spiking.(Vzoom), linewidth=1.5, size=(300, 200), color=:gray30, label="")
vline!([STG_Vth], linewidth=1.5, color=:firebrick1, linestyle=:dashdot, label="V_th")
p3 = plot(V, STG_gu_spiking.(V), linewidth=1.5, legend=false, size=(600, 200), color=:gray30)
vline!([Vzoom[1]], linewidth=1.5, color=:black, linestyle=:dash)
vline!([Vzoom[end]], linewidth=1.5, color=:black, linestyle=:dash)
xlabel!("V")
ylabel!("g_u")
p3zoom = plot(Vzoom, STG_gu_spiking.(Vzoom), linewidth=1.5, legend=false, size=(300, 200), color=:gray30)
vline!([STG_Vth], linewidth=1.5, color=:firebrick1, linestyle=:dashdot)
xlabel!("V")
CC = plot(p1, p1zoom, p2, p2zoom, p3, p3zoom, layout=(3, 2), size=(900, 600), margins=10px, dpi=500)
```
![](https://raw.githubusercontent.com/arthur-fyon/NmodController.jl/main/docs/src/assets/NmodController_DICs_STGspiking.png)

The next few lines of code show how to compute and plot DICs, sensitivity matrix and threshold voltage of the bursting STG model using *NmodController.jl* for the sake of comparison.

```julia
# Defining some timescales
STG_tauFast = STG_tau_mNa
STG_tauSlow = STG_tau_mKd
STG_tauUltraslow = STG_tau_mH

# Computing DICs for the bursting model
STG_gf_bursting, STG_gs_bursting, STG_gu_bursting = computeDICs(STG_bursting, STG_tauFast, 
    STG_tauSlow, STG_tauUltraslow, tauCa=500000.)

# Plotting
V = -80 : 0.1 : 40.
Vzoom = -53. : 0.01 : -48.
p1 = plot(V, STG_gf_bursting.(V), linewidth=1.5, legend=false, size=(600, 200), color=:gray30)
vline!([Vzoom[1]], linewidth=1.5, color=:black, linestyle=:dash)
vline!([Vzoom[end]], linewidth=1.5, color=:black, linestyle=:dash)
ylabel!("g_f")
p1zoom = plot(Vzoom, STG_gf_bursting.(Vzoom), linewidth=1.5, legend=false, size=(300, 200), color=:gray30)
vline!([STG_Vth], linewidth=1.5, color=:firebrick1, linestyle=:dashdot)
p2 = plot(V, STG_gs_bursting.(V), linewidth=1.5, legend=false, size=(600, 200), color=:gray30)
vline!([Vzoom[1]], linewidth=1.5, color=:black, linestyle=:dash)
vline!([Vzoom[end]], linewidth=1.5, color=:black, linestyle=:dash)
ylabel!("g_s")
p2zoom = plot(Vzoom, STG_gs_bursting.(Vzoom), linewidth=1.5, size=(300, 200), color=:gray30, label="")
vline!([STG_Vth], linewidth=1.5, color=:firebrick1, linestyle=:dashdot, label="V_th")
p3 = plot(V, STG_gu_bursting.(V), linewidth=1.5, legend=false, size=(600, 200), color=:gray30)
vline!([Vzoom[1]], linewidth=1.5, color=:black, linestyle=:dash)
vline!([Vzoom[end]], linewidth=1.5, color=:black, linestyle=:dash)
xlabel!("V")
ylabel!("g_u")
p3zoom = plot(Vzoom, STG_gu_bursting.(Vzoom), linewidth=1.5, legend=false, size=(300, 200), color=:gray30)
vline!([STG_Vth], linewidth=1.5, color=:firebrick1, linestyle=:dashdot)
xlabel!("V")
CC = plot(p1, p1zoom, p2, p2zoom, p3, p3zoom, layout=(3, 2), size=(900, 600), margins=10px, dpi=500)
```
![](https://raw.githubusercontent.com/arthur-fyon/NmodController.jl/main/docs/src/assets/NmodController_DICs_STGbursting.png)

# Computing DICs and threshold voltage on the DA model
Once the DA model is initialized, it is very straightforward to compute DICs, sensitivity matrix or threshold voltage. For comparison, two DA models with different maximum ion channel conductances will be used. `DA_spiking` that exhibits a tonic spiking behavior with
- sodium current: $\bar{g}_\mathrm{Na} =$ 30;
- delayed rectified potassium current: $\bar{g}_\mathrm{Kd} =$ 5;
- L-type calcium current: $\bar{g}_\mathrm{CaL} =$ 0.03;
- N-type calcium current: $\bar{g}_\mathrm{CaN} =$ 0.03;
- ERG current: $\bar{g}_\mathrm{ERG} =$ 0.12;
- NMDA current: $\bar{g}_\mathrm{NMDA} =$ 0;
- leakage current: $g_\mathrm{leak} =$ 0.01;

and `DA_bursting` that exhibits a bursting behavior with
- sodium current: $\bar{g}_\mathrm{Na} =$ 30;
- delayed rectified potassium current: $\bar{g}_\mathrm{Kd} =$ 5;
- L-type calcium current: $\bar{g}_\mathrm{CaL} =$ 0.02;
- N-type calcium current: $\bar{g}_\mathrm{CaN} =$ 0.15;
- ERG current: $\bar{g}_\mathrm{ERG} =$ 0.12;
- NMDA current: $\bar{g}_\mathrm{NMDA} =$ 0;
- leakage current: $g_\mathrm{leak} =$ 0.01.

Note that the NMDA current conductance is set to 0 since it is not involved in the intrinsic neuronal excitability, given that it involves specific NMDA receptors. 

The next few lines of code show how to compute and plot DICs, sensitivity matrix and threshold voltage of the spiking DA model using *NmodController.jl*.

```julia
# Defining some timescales
DA_tauFast = DA_tau_mNa
DA_tauSlow = DA_tau_mKd
DA_tauUltraslow(V) = 100.

# Computing DICs for the spiking model
DA_gf_spiking, DA_gs_spiking, DA_gu_spiking = computeDICs(DA_spiking, DA_tauFast, DA_tauSlow, DA_tauUltraslow, Mg=1.4)

# Computing the sensitivity matrix (which is a voltage function that can be called by DA_S_spiking(V))
DA_S_spiking = computeDICs(DA_spiking, DA_tauFast, 
    DA_tauSlow, DA_tauUltraslow, onlyS=true, Mg=1.4)

# Computing the threshold voltage
DA_Vth = -55. # Here, the algorithm is not able to find any threshold voltage

# Plotting
V = -80 : 0.1 : 40.
Vzoom = -58. : 0.01 : -53.
p1 = plot(V, DA_gf_spiking.(V), linewidth=1.5, legend=false, size=(600, 200), color=:gray30)
vline!([Vzoom[1]], linewidth=1.5, color=:black, linestyle=:dash)
vline!([Vzoom[end]], linewidth=1.5, color=:black, linestyle=:dash)
ylabel!("g_f")
p1zoom = plot(Vzoom, DA_gf_spiking.(Vzoom), linewidth=1.5, legend=false, size=(300, 200), color=:gray30)
vline!([DA_Vth], linewidth=1.5, color=:firebrick1, linestyle=:dashdot)
p2 = plot(V, DA_gs_spiking.(V), linewidth=1.5, legend=false, size=(600, 200), color=:gray30)
vline!([Vzoom[1]], linewidth=1.5, color=:black, linestyle=:dash)
vline!([Vzoom[end]], linewidth=1.5, color=:black, linestyle=:dash)
ylabel!("g_s")
p2zoom = plot(Vzoom, DA_gs_spiking.(Vzoom), linewidth=1.5, size=(300, 200), color=:gray30, label="")
vline!([DA_Vth], linewidth=1.5, color=:firebrick1, linestyle=:dashdot, label="V_th")
p3 = plot(V, DA_gu_spiking.(V), linewidth=1.5, legend=false, size=(600, 200), color=:gray30)
vline!([Vzoom[1]], linewidth=1.5, color=:black, linestyle=:dash)
vline!([Vzoom[end]], linewidth=1.5, color=:black, linestyle=:dash)
xlabel!("V")
ylabel!("g_u")
p3zoom = plot(Vzoom, DA_gu_spiking.(Vzoom), linewidth=1.5, legend=false, size=(300, 200), color=:gray30)
vline!([DA_Vth], linewidth=1.5, color=:firebrick1, linestyle=:dashdot)
xlabel!("V")
CC = plot(p1, p1zoom, p2, p2zoom, p3, p3zoom, layout=(3, 2), size=(900, 600), margins=10px, dpi=500)
```
![](https://raw.githubusercontent.com/arthur-fyon/NmodController.jl/main/docs/src/assets/NmodController_DICs_DAspiking.png)

The next few lines of code show how to compute and plot DICs, sensitivity matrix and threshold voltage of the bursting DA model using *NmodController.jl* for the sake of comparison.

```julia
# Defining some timescales
DA_tauFast = DA_tau_mNa
DA_tauSlow = DA_tau_mKd
DA_tauUltraslow(V) = 100.

# Computing DICs for the bursting model
DA_gf_bursting, DA_gs_bursting, DA_gu_bursting = computeDICs(DA_bursting, DA_tauFast, DA_tauSlow, DA_tauUltraslow, Mg=1.4)

# Plotting
V = -80 : 0.1 : 40.
Vzoom = -58. : 0.01 : -53.
p1 = plot(V, DA_gf_bursting.(V), linewidth=1.5, legend=false, size=(600, 200), color=:gray30)
vline!([Vzoom[1]], linewidth=1.5, color=:black, linestyle=:dash)
vline!([Vzoom[end]], linewidth=1.5, color=:black, linestyle=:dash)
ylabel!("g_f")
p1zoom = plot(Vzoom, DA_gf_bursting.(Vzoom), linewidth=1.5, legend=false, size=(300, 200), color=:gray30)
vline!([DA_Vth], linewidth=1.5, color=:firebrick1, linestyle=:dashdot)
p2 = plot(V, DA_gs_bursting.(V), linewidth=1.5, legend=false, size=(600, 200), color=:gray30)
vline!([Vzoom[1]], linewidth=1.5, color=:black, linestyle=:dash)
vline!([Vzoom[end]], linewidth=1.5, color=:black, linestyle=:dash)
ylabel!("g_s")
p2zoom = plot(Vzoom, DA_gs_bursting.(Vzoom), linewidth=1.5, size=(300, 200), color=:gray30, label="")
vline!([DA_Vth], linewidth=1.5, color=:firebrick1, linestyle=:dashdot, label="V_th")
p3 = plot(V, DA_gu_bursting.(V), linewidth=1.5, legend=false, size=(600, 200), color=:gray30)
vline!([Vzoom[1]], linewidth=1.5, color=:black, linestyle=:dash)
vline!([Vzoom[end]], linewidth=1.5, color=:black, linestyle=:dash)
xlabel!("V")
ylabel!("g_u")
p3zoom = plot(Vzoom, DA_gu_bursting.(Vzoom), linewidth=1.5, legend=false, size=(300, 200), color=:gray30)
vline!([DA_Vth], linewidth=1.5, color=:firebrick1, linestyle=:dashdot)
xlabel!("V")
CC = plot(p1, p1zoom, p2, p2zoom, p3, p3zoom, layout=(3, 2), size=(900, 600), margins=10px, dpi=500)
```
![](https://raw.githubusercontent.com/arthur-fyon/NmodController.jl/main/docs/src/assets/NmodController_DICs_DAbursting.png)
