% -------------------------------------------------------------------------
% Parameter_identification_main.m
% -------------------------------------------------------------------------
% Parameter_identification_main identifies parameters using experimentally 
% measured current and voltage data.

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
% The input for PSO identification is [I, V, t, SOC].
% -------------------------------------------------------------------------
% In this toolbox, the NMC experiment data are
% stored in a folder named "0_Data".
% Each data file contains 5 types of data:
% (1) C/2 charge, (2) C/2 discharge, (3) C/40 charge, (4) C/40 discharge, (5) HPPC.
% Each file contains (1) Voltage, (2) Current, (3) Time from the diagnostic cycle.


tot_nmc_cell=13; % This means data from cell #13 will be used for identification.

for tot_nmc=1:1:length(tot_nmc_cell) 
clearvars -except tot_nmc tot_nmc_cell


% -------------------------------------------------------------------------
% Step 1: Identify stoichiometry values using the half-cell model.
% -------------------------------------------------------------------------
% Parameter Identification Based on Up and Un Curve Fitting
% The basic assumption is that under very low C-rates, battery behavior follows thermodynamic principles. 
% Therefore, V = Up(θp) - Un(θn).
% Using this relationship, the stoichiometry parameters [θp, θn] can be identified.
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Because the cell used in the NMC experiment is Si+Gr, we have hysteresis.
% Therefore, we use Up_cc and Un_cc for charge data
% and use Up_dis and Un_dis for discharge data.
% This means we need to identify 4 stoichiometry [θp, θn] for Up_cc and Un_cc
% and identify 4 stoichiometry [θp, θn] for Up_dis and Un_dis.
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Load the Up and Un curves. 
% They represent the relationship between the state of lithium (θ) and the open-circuit voltage (OCV).
% The x-axis for both Up and Un is the state of lithium (θ).
% The y-axis for both Up and Un is the open-circuit voltage (OCV).
% -------------------------------------------------------------------------
load('OCVdatare.mat')
data.OCV=OCV;

% -------------------------------------------------------------------------
% Prepare the [I, V, t, SOC] data.
% C/40 charge and discharge will be used
% -------------------------------------------------------------------------
cell_num=tot_nmc_cell(tot_nmc);
if cell_num<10
FullVraw=load(['../0_Data/NMC-00',num2str(cell_num),'.mat']);
else
FullVraw=load(['../0_Data/NMC-0',num2str(cell_num),'.mat']);
end

% -------------------------------------------------------------------------
% The sampling frequency of the raw experimental data might be different.
% Here, we add a data preprocessing step to resample it to 1 Hz.
% -------------------------------------------------------------------------
% Re-sampling C/40 discharge data
int_num_DC=floor(FullVraw.co40dischargetime(end));
int_time_DC=1:int_num_DC;
int_V_DC=interp1(FullVraw.co40dischargetime',FullVraw.co40dischargeV',int_time_DC);
int_I_DC=interp1(FullVraw.co40dischargetime',FullVraw.co40dischargeI',int_time_DC);
% Re-sampling C/40 charge data
int_num_CC=floor(FullVraw.co40chargetime(end));
int_time_CC=1:int_num_CC;
int_V_CC=interp1(FullVraw.co40chargetime',FullVraw.co40chargeV',int_time_CC);
int_I_CC=interp1(FullVraw.co40chargetime',FullVraw.co40chargeI',int_time_CC);

% -------------------------------------------------------------------------
% Using coulomb counting to calculate SOC for C/40 discharge data.
cap_DC=abs(sum(int_I_DC)/3600);
Ifull_DC=int_I_DC;
socfull_DC(1)=1; 
for kk=1:1:size(Ifull_DC,2)-1
socfull_DC(kk+1)=socfull_DC(kk)+Ifull_DC(kk)/3600/cap_DC;
end
socfull_DC=socfull_DC'*100;
Vocfull_DC=int_V_DC;
num_DC=size(Vocfull_DC,2); 
data.Vocfull_DC=Vocfull_DC;
data.socfull_DC=socfull_DC;
data.intnum_DC=num_DC;

% -------------------------------------------------------------------------
% Using coulomb counting to calculate SOC for C/40 charge data.
cap_CC=abs(sum(int_I_CC)/3600);
Ifull_CC=int_I_CC;
socfull_CC(1)=0;
for kk=1:1:size(Ifull_CC,2)-1
socfull_CC(kk+1)=socfull_CC(kk)+Ifull_CC(kk)/3600/cap_CC;
end
socfull_CC=socfull_CC'*100;
Vocfull_CC=int_V_CC;
socfull_CC=flip(socfull_CC);
Vocfull_CC=flip(Vocfull_CC);
num_CC=size(Vocfull_CC,2); 
data.Vocfull_CC=Vocfull_CC;
data.socfull_CC=socfull_CC;
data.intnum_CC=num_CC;
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Particle swarm optimization (PSO) is used for identification
% This part is the PSO settings
% -------------------------------------------------------------------------
% Setting Upper and Lower bounds for 8 θ values
% Parameters 1-4 are for discharge: [θp100, θp0, θn100, θn0]
% Parameters 5-8 are for charge: [θp100, θp0, θn100, θn0] 
% -------------------------------------------------------------------------
Initial_position=[0.049,0.84,0.98,0.032,0.075,0.85,0.95,0.007]; % Initial guess values
Upper_bound=[0.4 1 1 0.4 0.4 1 1 0.4]; % Set upper bound
Lower_bound=[0 0 0 0 0 0 0 0];         % Set lower bound

% -------------------------------------------------------------------------
% PSO algorithm settings
% -------------------------------------------------------------------------
swarm_size=20;   % Number of particles in PSO
max_stalls=200;  % If the best values don't change in 200 iterations, stop.

% -------------------------------------------------------------------------
% Run PSO to identify [θp100, θp0, θn100, θn0] and [θp100, θp0, θn100, θn0]
% -------------------------------------------------------------------------
% After identification, x_opt contains identified 8 θ values
[x_opt, fval, history] = pso_for_step1(swarm_size, max_stalls, Initial_position, ...
                                       Lower_bound, Upper_bound,data);

% Parameters 1-4 are for discharge: [θp100, θp0, θn100, θn0]
thetap100=x_opt(1);
thetap0=x_opt(2);
thetan100=x_opt(3);
thetan0=x_opt(4);

% Parameters 5-8 are for charge: [θp100, θp0, θn100, θn0]
thetap100ch=x_opt(5);
thetap0ch=x_opt(6);
thetan100ch=x_opt(7);
thetan0ch=x_opt(8);

% -------------------------------------------------------------------------
% Save identified parameters. These values will be used by ESPM.
% -------------------------------------------------------------------------
stovar = [thetan100ch, thetap100ch, thetan0ch, thetap0ch];
save('sto_cc.mat', 'stovar');
save(['sto_cc', num2str(cell_num), '.mat'], 'stovar');
stovar = [thetan100, thetap100, thetan0, thetap0];
save('sto_dis.mat', 'stovar');
save(['sto_dis', num2str(cell_num), '.mat'], 'stovar');


% -------------------------------------------------------------------------
% Validation
% After identification, the identified parameters will be used for
% simulation. The simulated voltage will then be compared with the
% experimentally measured voltage for validation purposes.
% -------------------------------------------------------------------------

% Model-simulated voltage for C/40 charge
intsocp=linspace(thetap100ch,thetap0ch,data.intnum_CC);
intOCVp=interp1(data.OCV.pe_socch,data.OCV.pe_Upch,intsocp,'linear','extrap');
intsocn=linspace(thetan0ch,thetan100ch,data.intnum_CC);
intOCVn=interp1(data.OCV.ne_socch,data.OCV.ne_Unch,intsocn,'linear','extrap');
Vsim_cc=intOCVp-flip(intOCVn);
soc_cc=linspace(0,1,length(Vsim_cc));

% Model-simulated voltage for C/40 discharge
intsocp=linspace(thetap100,thetap0,data.intnum_DC);
intOCVp=interp1(data.OCV.pe_socdis,data.OCV.pe_Updis,intsocp,'linear','extrap');
intsocn=linspace(thetan0,thetan100,data.intnum_DC);
intOCVn=interp1(data.OCV.ne_socdis,data.OCV.ne_Undis,intsocn,'linear','extrap');
Vsim_dis=intOCVp-flip(intOCVn);
soc_dis=linspace(0,1,length(Vsim_dis));

% -------------------------------------------------------------------------
% Reference values from experiment
% -------------------------------------------------------------------------
Vexp1=data.Vocfull_DC'; % Discharge data
Vexp2=data.Vocfull_CC'; % Charge data

% -------------------------------------------------------------------------
% Plot to compare model simulation and experimental data
% -------------------------------------------------------------------------
figure;
subplot(1,2,1)
plot(soc_cc,Vsim_cc,'Linewidth',2)
hold on
plot(soc_cc,Vexp2,'-.','Linewidth',2)
legend('Sim','Exp')
xlabel('SOC')
ylabel('OCP')
set(gca,'linewidth',1,'fontsize',14,'fontname','Arial');
subplot(1,2,2)
plot(soc_dis,Vsim_dis,'Linewidth',2)
hold on
plot(soc_dis,Vexp1,'-.','Linewidth',2)
legend('Sim','Exp')
xlabel('SOC')
ylabel('OCP')
set(gca,'linewidth',1,'fontsize',14,'fontname','Arial');
% Save figure as .fig file
savefig(gcf, '../Figures/Half_cell.fig'); 


% -------------------------------------------------------------------------
% Step 2: Identify the rest of the parameter values using the ESPM model.
% -------------------------------------------------------------------------
clearvars -except cell_num tot_nmc tot_nmc_cell
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
% Based on sensitivity and correlation analysis results (LSA_CA_main.m).
% C/40 discharge, C/2 discharge and HPPC data will be used
% -------------------------------------------------------------------------
tot_profile_flag=[0,2,3];

% -------------------------------------------------------------------------
% Set initial values for PSO
% pso_flag == 1: This is the first time to run PSO, and we use default parameter values as the initial values.
% pso_flag == 2: This is a subsequent run of PSO, and we use previous PSO parameter values as the initial values.
% -------------------------------------------------------------------------
% Use default values to run PSO for the first time. 
% In subsequent runs, the identified values from the previous run will be used as the initial values.
tot_pso_flag=[1,2,2]; 

tot_err=[];
for tot_use=1:1:3
clearvars -except  tot_profile_flag tot_pso_flag tot_use cell_num tot_err tot_nmc tot_nmc_cell
% -------------------------------------------------------------------------
% Load the raw data
% -------------------------------------------------------------------------
profile_flag=tot_profile_flag(tot_use);
if cell_num<10
data=load(['../0_Data/NMC-00',num2str(cell_num),'.mat']);
else
data=load(['../0_Data/NMC-0',num2str(cell_num),'.mat']);
end

if profile_flag==0
rawI=data.co40dischargeI;
rawt=data.co40dischargetime;
rawV=data.co40dischargeV;
elseif profile_flag==1
rawI=data.co2chargeI;
rawt=data.co2chargetime;
rawV=data.co2chargeV;
elseif profile_flag==2
rawI=data.co2dischargeI;
rawt=data.co2dischargetime;
rawV=data.co2dischargeV;
elseif profile_flag==3
rawI=data.hppc_I;
dele_index=max(find(rawI<-2.25 & rawI>-2.7));
rawI=data.hppc_I(1:dele_index);
rawt=data.hppc_time(1:dele_index);
rawV=data.hppc_V(1:dele_index);
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
cap_t = data.co40dischargetime;
cap_I = data.co40dischargeI;
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

% Set initial values for PSO
% pso_flag == 1: This is the first time to run PSO, and we use default parameter values as the initial values.
% pso_flag == 2: This is a subsequent run of PSO, and we use previous PSO parameter values as the initial values.
pso_flag=tot_pso_flag(tot_use);

if pso_flag==1
load(['re_iden.mat'])        % This is the defalut values

if profile_flag==0
    sto=load('sto_dis.mat');
    sto=sto.stovar;
    x_opt(5:8)=sto;
elseif profile_flag==1
    sto=load('sto_cc.mat');
    sto=sto.stovar;
    x_opt(5:8)=sto;   
elseif profile_flag==2
    sto=load('sto_dis.mat');
    sto=sto.stovar;
    x_opt(5:8)=sto;
elseif profile_flag==3
    sto=load('sto_dis.mat');
    sto=sto.stovar;
    x_opt(5:8)=sto;
end

elseif pso_flag==2   % This is the previous PSO values
load(['check.mat'])

if profile_flag==0
    sto=load('sto_dis.mat');
    sto=sto.stovar;
    x_opt(5:8)=sto;
elseif profile_flag==1
    sto=load('sto_cc.mat');
    sto=sto.stovar;
    x_opt(5:8)=sto;   
elseif profile_flag==2
    sto=load('sto_dis.mat');
    sto=sto.stovar;
    x_opt(5:8)=sto;
elseif profile_flag==3
    sto=load('sto_dis.mat');
    sto=sto.stovar;
    x_opt(5:8)=sto;
end

end

I_data=refall.I;
t_data=refall.t_data;

% -------------------------------------------------------------------------
% Set initial SOC values for different types of data.
% profile_flag = 0;   % C/40 discharge
% profile_flag = 1;   % C/2 charge
% profile_flag = 2;   % C/2 discharge
% profile_flag = 3;   % HPPC discharge
% -------------------------------------------------------------------------

if profile_flag==0
    SOC_IC=1; % C/40 discharge
elseif profile_flag==1
    SOC_IC=0; % C/2 charge
elseif profile_flag==2
    SOC_IC=1; % C/2 discharge
elseif profile_flag==3
    SOC_IC=1; % HPPC discharge
end

% -------------------------------------------------------------------------
% Prepare data for PSO
% -------------------------------------------------------------------------
data.model_type=model_type;
data.int_method=int_method;
data.I_data=I_data;
data.t_data=t_data;
data.SOC_IC=SOC_IC;
data.refV=refall.V;
data.refsocp=refall.soc_p;
data.refsocn=refall.soc_n;
data.profile_flag=profile_flag;
data.x_opt=x_opt;

% -------------------------------------------------------------------------
% In this toolbox, we can identify the following parameters in ESPM
% Parameter 1-2: Rp, Rn - Particle radius
% Parameter 3-4: epsp, epsn - Volume fraction for active materials
% Parameter 5-6: Dsp, Dsn - Solid-phase diffusion coefficient
% Parameter 7-8: kp, kn - Reaction rate constant
% Parameter 9-10: csmaxp, csmaxn - Maximum lithium concentration
% Parameter 11-12: Lp, Ln - Electrode thickness
% Parameter 13: De - Electrolyte diffusion coefficient
% Parameter 14: ke - Electrolyte conductivity
% Parameter 15: A - Cell area
% Parameter 16: Rs - Lumped resistance
% Parameter 17-19: Brugg for positive, separator, and negative electrode
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Step1: Setting initial guess values for PSO
% -------------------------------------------------------------------------
if pso_flag==1     % Using default values as initial guess values 
Initial_position=[5.7,5.9,0.71,0.76,15.0,15.2,2.2,5.95,4.52,4.63,-4.07,-4.12,-10.12,1.18,0.15,-1.75,1.69,1.0,1.06];
elseif pso_flag==2 % Using previous step identification results as initial guess values 
Initial_position=load('pre_pso.mat');
Initial_position=Initial_position.x_opt_iden;
end

% -------------------------------------------------------------------------
% Step2: Setting upper and lower bounds for identified parameters.
% Here, we first set the general bounds for all 19 parameters.
% -------------------------------------------------------------------------

if pso_flag==1 % This means it's the first time to run PSO, so we set the upper and lower bounds according to experience.
Upper_bound=[10,10,0.9,0.9,16,16,6,10,6,6,-2,-2,-5,2,0.3,-1,4,4,4];
Lower_bound=[4.4,4.4,0.5,0.5,12,12,2,2,2,2,-6,-6,-11,0.1,0.01,-3,1,1,1];

elseif pso_flag==2 % We set the upper and lower bounds according to previous identification results.
                   % More specifically, the lower bound is 10 times smaller
                   % and the upper bound is 10 times larger.

   Upper_bound= Initial_position;
   Lower_bound= Initial_position;
   
   for kk=1:1:size(Initial_position,2)
    if (Upper_bound(kk)>0)
       Upper_bound(kk)=Upper_bound(kk)*10;
    else
        Upper_bound(kk)=Upper_bound(kk)*0.1;
    end
    
    if (Lower_bound(kk)>0)
       Lower_bound(kk)=Lower_bound(kk)*0.1;
    else
        Lower_bound(kk)=Lower_bound(kk)*10;
    end
   end

end

% -------------------------------------------------------------------------
% If we want to fix a specific parameter, 
% the most straightforward way is to set the lower bound equal to the upper bound.
% -------------------------------------------------------------------------
% 1. Fix csmaxp and csmaxn
Upper_bound(9:10)=[log10(29583),log10(51760)];
Lower_bound(9:10)=[log10(29583),log10(51760)];

% 2. Fix Lp and Ln
Upper_bound(11:12)=[log10(85.2e-6),log10(75.6e-6)];
Lower_bound(11:12)=[log10(85.2e-6),log10(75.6e-6)];

% 3. Fix Bruggp, Bruggs, and Bruggn
Upper_bound(17:19)=[1.5,1.5,1.5];
Lower_bound(17:19)=[1.5,1.5,1.5];

% -------------------------------------------------------------------------
% Additionally, if we want to adjust the bounds for a specific parameter, we can redefine it.
% -------------------------------------------------------------------------
% 4. Set a very small range for the cell area
Upper_bound(15)=[0.16];
Lower_bound(15)=[0.15];

% 5. Set parameter range for electrolyte diffusivity De
Upper_bound(13)=[-10];
Lower_bound(13)=[-14];

% 6. Set parameter range for lumped resistance Rs
Upper_bound(16)=[log10(0.04)];
Lower_bound(16)=[log10(0.015)];

% -------------------------------------------------------------------------
% 7. Design the identification strategy for Rp and Rn
% Here, the following settings mean Rp and Rn are identified using
% profile 2 (C/2 discharge). Then, when we identify other parameters using
% profile 3 (HPPC), Rp and Rn values are fixed.
% -------------------------------------------------------------------------
if profile_flag ==0 || profile_flag ==1|| profile_flag ==2
Upper_bound(1:2)=[-log10(1e-6),-log10(1e-6)];
Lower_bound(1:2)=[-log10(8e-6),-log10(8e-6)];
else
Upper_bound(1:2)=[Initial_position(1),Initial_position(2)];
Lower_bound(1:2)=[Initial_position(1),Initial_position(2)];    
end

% -------------------------------------------------------------------------
% 8. Design the identification strategy for eps_p and eps_n
% Here, the following settings mean eps_p and eps_n are identified using
% profile 2 (C/2 discharge). Then, when we identify other parameters using
% profile 3 (HPPC), eps_p and eps_n values are fixed.
% -------------------------------------------------------------------------
if profile_flag ==0 || profile_flag ==1|| profile_flag ==2
Upper_bound(3:4)=[0.9,0.9];
Lower_bound(3:4)=[0.3,0.3];
else
Upper_bound(3:4)=[Initial_position(3),Initial_position(4)];
Lower_bound(3:4)=[Initial_position(3),Initial_position(4)];
end


% -------------------------------------------------------------------------
% PSO algorithm settings
% -------------------------------------------------------------------------
swarm_size=10;   % Number of particles in PSO
max_stalls=200;  % If the best values don't change in 200 iterations, stop.

% -------------------------------------------------------------------------
% Run PSO to identify ESPM model parameters
% -------------------------------------------------------------------------
% After identification, x_opt contains identified parameters
[x_opt_iden, fval, history] = pso_for_step2(swarm_size, max_stalls, Initial_position, ...
    Lower_bound, Upper_bound,data);

% After each step of identification, store the error and parameter values.
tot_err(tot_use)=fval;
save('pre_pso.mat','x_opt_iden')


% -------------------------------------------------------------------------
% Validation
% After identification, the identified parameters will be used for
% simulation. The simulated voltage and SOC will then be compared with the
% experimentally measured voltage for validation purposes.
% -------------------------------------------------------------------------

% Load identified parameters
x_opt(1:2)=x_opt_iden(1:2);
x_opt(3:4)=x_opt_iden(3:4);
x_opt(9:10)=x_opt_iden(5:6);
x_opt(11:12)=x_opt_iden(7:8);
x_opt(15:16)=x_opt_iden(9:10);
x_opt(17:18)=x_opt_iden(11:12);
x_opt(20)=x_opt_iden(13);
x_opt(22)=x_opt_iden(14);
x_opt(23)=x_opt_iden(15);
x_opt(28)=x_opt_iden(16);
x_opt(25:27)=x_opt_iden(17:19);

% Model Simulation
T_amb=23;
[V_cell, R_l, T_core, T_surf,soc_bulk_n, soc_bulk_p, cs_n, cs_p,...
          V_oc,param,ocp_p,ocp_n,eta_p,eta_n,ce,eta_ele] = ESPM_main(x_opt,t_data,I_data,SOC_IC,T_amb,profile_flag);

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
savefig(gcf, ['../Figures/Validation_',num2str(tot_use),'.fig']); % Save figure as .fig file

save('check.mat','x_opt')
end

% -------------------------------------------------------------------------
% Save identified parameters
% -------------------------------------------------------------------------

save(['iden_par',num2str(cell_num),'.mat'],'x_opt_iden')
save(['iden_par_err',num2str(cell_num),'.mat'],'tot_err')

end

% -------------------------------------------------------------------------
% Show all the results of identification
% Fig.1: Fitting results using the Half-cell model
% Fig.2: Fitting results using ESPM under profile 0 (C/40 discharge)
% Fig.3: Fitting results using ESPM under profile 2 (C/2 discharge)
% Fig.4: Fitting results using ESPM under profile 3 (HPPC)
% -------------------------------------------------------------------------
close all
clear all
clc

open ('../Figures/Half_cell.fig')
open ('../Figures/Validation_1.fig')
open ('../Figures/Validation_2.fig')
open ('../Figures/Validation_3.fig')
