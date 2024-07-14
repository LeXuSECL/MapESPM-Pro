% -------------------------------------------------------------------------
% Battery_Simulator.m
% -------------------------------------------------------------------------
% Battery_Simulator is used for ESPM model simulation.
% After identifying ESPM model parameters (use Parameter_identification_main.m),
% the identified ESPM model parameters are stored in the "ESPM_Parameters" folder.
% By running this script, the simulated Voltage and SOC can be compared
% with the experimentally measured values.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MapESPM-Pro: Mass-Preserving Extended Single Particle Model with Parameter Optimization for Control-Oriented Applications

% Copyright (c) 2024 Mass-Preserving Extended Single Particle Model with Parameter Optimization for Control-Oriented Applications (MapESPM-Pro)
% MapESPM-Pro is freely distributed software under the MIT License 
% v1.0.0.: Released 07/2024 

% Written by Le Xu (lexu1209@stanford.edu)
% PI: Prof. Simona Onori (sonori@stanford.edu)
% Stanford Energy Control Group, Stanford University
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;clear;close all
addpath ('../2_Functions')

% -------------------------------------------------------------------------
% The input for Battery_Simulator.m are:
% (1) Time
% (2) Current
% (3) Experimental measured Voltage (For comparsion use)
% (4) Experimental measured cell SOC (For comparsion use)
% -------------------------------------------------------------------------
% In this toolbox, the ZOOX aging experiment data conducted at SLAC are
% stored in a folder named "0_Data".
% Each data file contains 5 types of data:
% (1) C/2 charge, (2) C/2 discharge, (3) C/40 charge, (4) C/40 discharge, (5) HPPC.
% Each file contains (1) Voltage, (2) Current, (3) Time from the diagnostic cycle.


tot_use_cell=1; % This means data from LG cell #1 will be used for simulation.

for current_use_cell=1:1:length(tot_use_cell) 
clearvars -except current_use_cell tot_use_cell

cell_num=tot_use_cell(current_use_cell);


clearvars -except cell_num current_use_cell tot_use_cell
addpath ('../2_Functions')
addpath('../1_Solver/casadi-3.6.5-windows64-matlab2018b')
import casadi.*

% -------------------------------------------------------------------------
% For each cell data file, we have the following types of data:
% profile_flag = 0;   % C/40 discharge
% profile_flag = 1;   % C/2 charge
% profile_flag = 2;   % C/2 discharge
% profile_flag = 3;   % HPPC discharge
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% The following current input will be used for this Battery Simulator.
% (1) C/40 discharge   (profile_flag = 0) 
% (2) C/2 discharge    (profile_flag = 2) 
% (3) HPPC data        (profile_flag = 3) 
% -------------------------------------------------------------------------
tot_profile_flag=[0,0,3];

for tot_use=1:1:size(tot_profile_flag,2)
clearvars -except  tot_profile_flag tot_pso_flag tot_use cell_num tot_err current_use_cell tot_use_cell
% -------------------------------------------------------------------------
% Load the raw data
% -------------------------------------------------------------------------
profile_flag=tot_profile_flag(tot_use);
if cell_num<10
data=load(['../0_Data/NMC_LG/LG-00',num2str(cell_num),'.mat']);
else
data=load(['../0_Data/NMC_LG/LG-0',num2str(cell_num),'.mat']);
end

if profile_flag==0
rawI=data.co20dischargeI;
rawt=data.co20dischargetime;
rawV=data.co20dischargeV;
elseif profile_flag==1

elseif profile_flag==2

elseif profile_flag==3
rawI=data.hppc_I;
rawt=data.hppc_time;
rawV=data.hppc_V;
end

% -------------------------------------------------------------------------
% The sampling frequency of the raw experimental data might be different.
% Here, we add a data preprocessing step to resample it to 1 Hz.
% -------------------------------------------------------------------------
intx=1:1:floor(rawt(end));
refall.I=interp1(rawt,rawI,intx);
refall.I=-refall.I;
refall.V=interp1(rawt,rawV,intx);
refall.t_data=intx;

% -------------------------------------------------------------------------
% Using the Coulomb counting method to calculate battery SOC.
% For Coulomb counting method, battery capacity is needed.
% Here, the battery capacity is obtained from the C/40 data.
% -------------------------------------------------------------------------
% Step1 Obtain battery capacity 
cap_t = data.co20dischargetime;
cap_I = data.co20dischargeI;
intx = 1:1:floor(cap_t(end));
cap_I = interp1(cap_t, cap_I, intx);
cap = 0;

for kk = 1:1:size(cap_I, 2)
  cap = cap + abs(cap_I(kk)) / 3600;
end

% Step2 Coulomb counting method
% -------------------------------------------------------------------------
% Set the initial SOC value for Coulomb counting
% -------------------------------------------------------------------------

if profile_flag == 0      % C/40 discharge
    soc(1) = 1;
elseif profile_flag == 1  % C/2 charge
    soc(1) = 0;
elseif profile_flag == 2  % C/2 discharge
    soc(1) = 1;
elseif profile_flag == 3  % HPPC discharge
    soc(1) = 1;
end

for kk = 2:size(refall.I, 2)
    soc(kk) = soc(kk - 1) - refall.I(kk) / 3600 / cap;  
end

refall.soc_p = soc;
refall.soc_n = soc;

% -------------------------------------------------------------------------
% ESPM model settings

% In this code, we use FVM to solve all the governing equations.
% Using FVM, we need to further calculate surface concentration (cssurf).
% In FVM, we only obtain average concentration at each control volume (CV).
% Here, we use different interpolation methods to calculate FVM cssurf.
% More details about the FVM scheme can be found here:
% Xu, L., Cooper, J., Allam, A., Onori, S., 
% "Comparative Analysis of Numerical Methods for Lithium-Ion Battery Electrochemical Modeling," 
% Journal of The Electrochemical Society, 170 (12), 120525, 2023.
% -------------------------------------------------------------------------
% Define interpolation method 
% FVM-S1: linear interpolation method to calculate surface concentration
% FVM-S2: Hermite interpolation method to calculate surface concentration
% int_method = 1; % linear
% int_method = 2; % Hermite

int_method = 1;

% Define using SPM or ESPM 
% model_type = 1; % SPM
% model_type = 2; % ESPM
model_type = 2;

% Define which type of battery is used
% battery_type=1; % NMC_Pansonic
% battery_type=2; % NMC_LG
battery_type=2;




I_data=refall.I;
t_data=refall.t_data;

% -------------------------------------------------------------------------
% Set initial SOC values for different types of data.
% profile_flag = 0;   % C/40 discharge
% profile_flag = 1;   % Nan
% profile_flag = 2;   % Nan
% profile_flag = 3;   % HPPC discharge
% -------------------------------------------------------------------------

if profile_flag==0
    SOC_IC=1; % C/40 discharge
elseif profile_flag==1

elseif profile_flag==2

elseif profile_flag==3
    SOC_IC=1; % HPPC discharge
end



% -------------------------------------------------------------------------
% Load ESPM model parameters
% -------------------------------------------------------------------------
x_opt=load('../ESPM_Parameters/NMC_LG_cell_1.mat');
x_opt=x_opt.x_opt;

% Model Simulation
T_amb=23;
[V_cell, R_l, T_core, T_surf,soc_bulk_n, soc_bulk_p, cs_n, cs_p,...
          V_oc,param,ocp_p,ocp_n,eta_p,eta_n,ce,eta_ele] = ESPM_main(x_opt,t_data,I_data,SOC_IC,T_amb,profile_flag,int_method,model_type,battery_type);

% Plot and Save figures      
figure;
subplot(2,2,[1,2])
plot(t_data,refall.V,'linewidth',2)
hold on
plot(t_data,V_cell,'-.','linewidth',2)
legend('Exp','Sim')
xlabel('Time [s]')
ylabel('Voltage [V]')
title(['RMSE=',num2str(rms(V_cell-refall.V)*1000),' mV'])
set(gca,'linewidth',1,'fontsize',14,'fontname','Arial');

subplot(2,2,3);
plot(t_data,refall.soc_p,'linewidth',2)
hold on
plot(t_data,soc_bulk_p,'-.','linewidth',2)
legend('Exp','Sim')
xlabel('Time [s]')
ylabel('SOC_p [-]')
ylim([0,1])
title(['RMSE=',num2str(rms(soc_bulk_p-refall.soc_p)*100),' %'])
set(gca,'linewidth',1,'fontsize',14,'fontname','Arial');

subplot(2,2,4);
plot(t_data,refall.soc_n,'linewidth',2)
hold on
plot(t_data,soc_bulk_n,'-.','linewidth',2)
legend('Exp','Sim')
xlabel('Time [s]')
ylabel('SOC_n [-]')
ylim([0,1])
title(['RMSE=',num2str(rms(soc_bulk_n-refall.soc_n)*100),' %'])
set(gca,'linewidth',1,'fontsize',14,'fontname','Arial');
set(gcf,'unit','centimeters','position',[5 5 26 20])
savefig(gcf, ['../Figures/Testing_',num2str(tot_use),'.fig']); % Save figure as .fig file

end



end

% -------------------------------------------------------------------------
% Show all the results of Simulation
% Fig.1: Fitting results using the Half-cell model
% Fig.2: Fitting results using ESPM under profile 0 (C/40 discharge)
% Fig.3: Fitting results using ESPM under profile 2 (C/2 discharge)
% Fig.4: Fitting results using ESPM under profile 3 (HPPC)
% -------------------------------------------------------------------------
close all
clear all
clc


open ('../Figures/Testing_1.fig')
open ('../Figures/Testing_2.fig')
open ('../Figures/Testing_3.fig')
