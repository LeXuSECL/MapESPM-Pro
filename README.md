# MapESPM-Pro: Mass-Preserving Extended Single Particle Model with Parameter Optimization for Control-Oriented Applications
This toolbox incorporates an identification algorithm that enables users to determine model parameters for any type of lithium-ion battery using current-voltage experimental data. By formulating an optimization problem, this toolbox minimizes the difference between simulated and experimentally measured voltage and state-of-charge (SOC) data to identify the parameters.

# Why MapESPM-Pro is needed?
The physics-based Electrochemical Single Particle Model (ESPM) offers significantly lower complexity compared to the Doyle-Fuller-Newman (DFN) model, making it suitable for on-board applications in Battery Management Systems (BMS). However, to accurately predict the dynamic behavior of actual batteries, the ESPM requires precise calibration of its parameters. Unfortunately, existing ESPM modeling tools lack the capability for parameter identification.
Commercially available COMSOL Multiphysics® software is widely used for simulating both the Doyle-Fuller-Newman (DFN) model and the Single Particle Model (SPM). By integrating COMSOL® with MATLAB® through LiveLink™, researchers can perform parameter identification for the DFN and SPM models. However, a similar framework for the Electrochemical Single Particle Model (ESPM) has yet to be proposed. Furthermore, the high licensing costs of COMSOL present a significant barrier to code accessibility and collaboration among researchers. On the other hand, several open-source ESPM model simulation tools, such as PyBaMM, SPMeT, and Spectral Li-ion SPM, are available. Despite their accessibility, none of these open-source tools currently offer parameter identification capabilities using experimental data.

In addition to parameter identification, mass-preserving properties are critical when implementing battery models for long-term real-world applications. A mass-conserved battery model ensures that there is no accumulated numerical error over time, which is essential for maintaining the accuracy and reliability of simulations. Among above mentioned software/tools, PyBaMM is the only software that solve ESPM model governing equations using a mass-preserving numerical finite volume method (FVM) scheme.

# What can MapESPM-Pro do?
- Local sensitivity analysis of model parameters
- Correlation analysis of model parameters
- Multi-step identification of model parameters based purely on current-voltage data
- Battery simulator for voltage and state of charge (SOC)

# Examples
You will find example codes that will help you get started.
- Sensitivity_and_Correlation_Analysis: examples showing how to perform sensitivity and correlation analysis for ESPM model.

  run "LSA_CA_main.m"
  
# Software dependencies
- MATLAB 2018b and later
- MATLAB Global Optimization Toolbox
- MATLAB Parallel Computing Toolbox
- CasADi
