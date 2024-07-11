
% -------------------------------------------------------------------------
% ModelParameters.m
% -------------------------------------------------------------------------
% ModelParameters is used to initialize ESPM model parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MapESPM-Pro: Mass-Preserving Extended Single Particle Model with Parameter Optimization for Control-Oriented Applications

% Copyright (c) 2024 Mass-Preserving Extended Single Particle Model with Parameter Optimization for Control-Oriented Applications (MapESPM-Pro)
% MapESPM-Pro is freely distributed software under the MIT License 
% v1.0.0.: Released 07/2024 

% Written by Le Xu (lexu1209@stanford.edu)
% PI: Prof. Simona Onori (sonori@stanford.edu)
% Stanford Energy Control Group, Stanford University
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Initialize Model Parameters (SI Units)
param.Rg = 8.314;                   % Gas constant
param.F = 96487;                    % Faraday's constant
param.alpha_cell = 0.5;             % Anode/Cathode transfer coefficient
param.ce0 = 2000;                   % Average electrolyte concentration [mol/m^3]
param.c_n_max = 10^(x(15));         % Maximum anode concentration [mol/m^3]
param.c_p_max = 10^(x(16));         % Maximum cathode concentration [mol/m^3]

%% GEOMETRIC PARAMETERS
param.Rs_n = 10^(-x(1));            % Anode particle radius
param.Rs_p = 10^(-x(2));            % Cathode particle radius         
param.A = x(23);                    % Area (averaged from anode/cathode area)
param.Ln = 10^(x(17));              % Anode thickness 
param.Lp = 10^(x(18));              % Cathode thickness 
param.epsilon_n = x(3);             % Anode solid phase volume fraction (1-poro)
param.epsilon_p = x(4);             % Cathode solid phase volume fraction (1-poro)
param.theta100_n = x(5);            % Anode stoichiometry at 100% SOC
param.theta100_p = x(6);            % Cathode stoichiometry at 100% SOC
param.theta0_n = x(7);              % Anode stoichiometry at 0% SOC
param.theta0_p = x(8);              % Cathode stoichiometry at 0% SOC

% Electrolyte porosity
param.eps_filler_n = x(13);         % Volume fraction for filler in negative electrode
param.eps_filler_p = x(14);         % Volume fraction for filler in positive electrode
param.eps_el_n_ref = 1-param.epsilon_n-param.eps_filler_n;
param.eps_el_p = 1-param.epsilon_p-param.eps_filler_p;
param.eps_el_s = x(21);                           
param.kep=x(22);                    %Electrolyte conductivity 
param.kes=x(22);                    %Electrolyte conductivity 
param.ken=x(22);                    %Electrolyte conductivity 
param.De=10^(x(20));                %Electrolyte diffusion cofficient 

% Specific interfacial (electroactive) surface area for anode and cathode
param.a_sn = 3*param.epsilon_n/param.Rs_n;   
param.a_sp = 3*param.epsilon_p/param.Rs_p; 


%% TRANSPORT/KINETIC PARAMETERS
% Diffusion coefficients
param.Dsn_ref = 10^(-x(9));              % Anode diffusion coefficient
param.Dsp_ref = 10^(-x(10));             % Cathode diffusion coefficient

% Reaction rates
param.kn_ref = 10^(-x(11))/param.F;      % Anode reaction rate constant 
param.kp_ref = 10^(-x(12))/param.F;      % Cathode reaction rate constant

% Lumped Contact Resistance
param.R_l = 10^x(28);                    % Lumped resistance

% Additional Electrolyte Parameters
param.t0 = x(24);                        % Transference Number 
param.Ls = 10^(x(19));                   % Separator Thickness [m] 

% Bruggeman coefficients 
param.brugg_n = x(25); 
param.brugg_s = x(26);
param.brugg_p = x(27);


%% THERMAL MODEL PARAMETERS
param.Tref = 298;                   % Reference Temperature [K]

