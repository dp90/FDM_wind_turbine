# FDM_wind_turbine
Application of the finite difference method to a wind turbine. 

The wind turbine is modeled as a Timoshenko beam embedded in soil en loaded by waves and wind. The rotor nacelle is modeled as a mass at the top end of the beam. A figure of the studied object is shown in the Powerpoint file "Schematic_representation_wind_turbine". 
The solution is found through application of the finite difference method and the Matlab ode45 solver applied to the state-space vector. The problem is thus solved in the time domain. 

The Matlab files perform the following functions:
1. Main.m contains input parameters, creates soil stiffness and damping matrices and plots the final results. (Takes about 20 minutes to run, so stresses in the turbine over time are shown in Stresses_turbine_over_time.png).
2. solve_statespace_vector.m solves at each time step the state-space vector of the system. 
3. extract_wind_speed.m extracts the wind speeds at each time step from Wind_signal.mat as input for the solve_statespace_vector.m.
4. calculate_wave_loads_JONSWAP.m calculates the water displacements, velocities and resulting wave loads from the JONSWAP wave spectrum. 
