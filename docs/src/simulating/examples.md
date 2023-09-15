# Example on existing models in the litterature

## Simulating the uncontrolled STG model
Once the STG model is initialized, simulated it can be simplified using *NmodController.jl*. At first, we will simulate the model `STG_spiking` that exhibits a tonic spiking behavior with
- sodium current: $\bar{g}_\mathrm{Na} =$ 4000;
- T-type calcium current: $\bar{g}_\mathrm{CaT} =$ 3;
- slow calcium current: $\bar{g}_\mathrm{CaS} =$ 4;
- A-type potassium current: $\bar{g}_\mathrm{A} =$ 175;
- calcium controlled potassium current: $\bar{g}_\mathrm{KCa} =$ 110;
- delayed rectified potassium current: $\bar{g}_\mathrm{Kd} =$ 137;
- H type current: $\bar{g}_\mathrm{H} =$ 0.3;
- leakage current: $g_\mathrm{leak} =$ 0.01.

The next few lines of code show how to simulate the spiking STG model using *NmodController.jl*.

