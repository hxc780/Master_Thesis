#### Script: FitModel.rmd #####
- Fits the model, using warp.alg(.), on either synthetic data (Signal_Synthetic.csv) or real data (Signal_Real.csv)
- Plots the fitted signal, fitted xi, fitted g(x), fitted envelope
- USER INPUT: 
	(Chunk 1): Select real signal or synthetic signal
	(Chunk 2): Optional - Select runtime parameters: n.rep (warp iterations), n.particles (SMC particles), n.mem (when to start "squeezing")

#### Script: Generate Synthetic Data.rmd ####
- Generate and saves (.csv) a synthetic signal simulated from the model
- Plots the simulated signal
- USER INPUT:
	(Chunk 2): Select model parameters

#### Other scripts ####
- Initialization.R: Implements Simulated Annealing and AMPD (peak detection) for initialization of parameters.
- Filter - SMC.R: Implement the warping algorithm using SMC as filter for \xi, including SMC-ensemble for initialization of psi, omega (i.e. multiple trajectories for different configurations of psi,omega)
- Estimation.R: Main-script which selects the relevant initialization procedure and filter (here only SMC), before fitting the model.  
- HelperFunctions.R: All "minor functions" used withn the scipts above, including simulation of a signal and simulation of xi.


Contact:
Lars Reiter Nielsen
lrn@math.ku.dk
lars.reiter.nielsen@gmail.com