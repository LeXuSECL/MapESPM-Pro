% -------------------------------------------------------------------------
% LSA_CA_main.m
% -------------------------------------------------------------------------
% LSA_CA_main performs local sensitivity analysis (LSA) and corrolation analysis (CA)
% for the ESPM model parameters. This script will run LSA and CA under
% given current input and output identifiable ESPM model parameters. 
% These identifiable ESPM model parameters are sensitive and are uncorrelated with each other.
% Therefore, these parameters can be identified based on experimental data using "XX" script.
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
%% addpath for model functions and numerical solver Casadi
addpath ('../2_Functions')
addpath('../1_Solver/casadi-3.6.5-windows64-matlab2018b')
import casadi.*

%--------------------------------------------------------------------------
% ESPM mdoel settings
%--------------------------------------------------------------------------
global useNr int_method model_type
% In this toolbox, the finite difference method (FVM) is used to solve the
% ESPM model euqations. % Here we need to define the number of Control Volumes (CV)
useNr=21;

% In ESPM, we need solid-phase surface concentration (cssurf) to calculate battery voltage.
% In FVM, we only obtain average concentration at each CV.
% Here, we use two different interpolation methods to calculate FVM cssurf.
% ========== Define interpolation method ==========
% FVM-S1: linear interpolation method to calculate surface concentration
% FVM-S2: Hermite interpolation method to calculate surface concentration
% int_method = 1; % linear
% int_method = 2; % Hermite
% More details about the FVM scheme can be found here:
% Xu, L., Cooper, J., Allam, A., Onori, S., "Comparative Analysis of Numerical Methods for Lithium-Ion Battery Electrochemical Modeling," Journal of The Electrochemical Society, 170 (12), 120525, 2023.
int_method=1;

% =========Define using SPM or ESPM
% model_type=1; % SPM
% model_type=2; % ESPM
model_type=2;

%--------------------------------------------------------------------------
% Define data used for ESPM parameter sensitivty and corrolation analysis
%--------------------------------------------------------------------------
data=load('LSA_CA_Data.mat'); % All avaliable data are stored here
% For each cell data file, we have the following types of data:
% profile_flag = 0;   % C/40 discharge
% profile_flag = 1;   % C/2 charge
% profile_flag = 2;   % C/2 discharge
% profile_flag = 3;   % HPPC discharge
tot_profile_flag=[0,2,3]; % C/40 discharge, C/2 discharge and HPPC data will be used for analysis
SA_flag=1;
SA_results=[];
CA_results=[];
for totkk=1:1:3
profile_flag=tot_profile_flag(totkk);
% Load data
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

%--------------------------------------------------------------------------
% NOTE: The sampling frequency of the raw experimental data might be different.
% Here, we add a data preprocessing step to resample it to 1 Hz.
%--------------------------------------------------------------------------
intx=1:1:floor(rawt(end));
refall.I=interp1(rawt,rawI,intx);
refall.I=-refall.I;
refall.V=interp1(rawt,rawV,intx);
refall.t_data=intx;

%--------------------------------------------------------------------------
% Set Initial SOC [-]
%--------------------------------------------------------------------------
if profile_flag==0
    SOC_IC=1;            % C/40 discharge
elseif profile_flag==1
    SOC_IC=0;            % C/2 charge
elseif profile_flag==2
    SOC_IC=1;            % C/2 discharge
elseif profile_flag==3
    SOC_IC=1;            % HPPC discharge
end

%--------------------------------------------------------------------------
% Prepare data for LSA and CA 
%--------------------------------------------------------------------------
data.I_data=refall.I;
data.t_data=refall.t_data;
data.SOC_IC=SOC_IC;


%--------------------------------------------------------------------------
% In total, we analyze the sensitivity of 13 parameters.
% We use the Fisher information matrix to calculate sensitivity.
% The first step is to generate values for the sensitivity matrix by perturbing each parameter one at a time
% Details can be found in
% https://iopscience.iop.org/article/10.1149/1945-7111/ad1293 (Eq. 50)
%--------------------------------------------------------------------------
% Perturbation coefficient for LSA [-]
SA_delta=0.1;  % Set 10% perturbation

%--------------------------------------------------------------------------
% 1. SA_matrix is the calculated sensitivity matrix.
% 2. index_SA is the sensitivity ranking for the 13 ESPM parameters in descending order.
% 3. gsax_labels are the labels of the 13 ESPM parameters corresponding to the sensitivity ranking.
% 4. SA_ESPM contains the sensitivity values of the ESPM parameters.
%--------------------------------------------------------------------------
[x,y,gsax_labels,SA_matrix,index_SA,SA_ESPM]=function_LSA_CA(SA_delta,data,profile_flag);
SA_results{totkk}=[index_SA;y];


%--------------------------------------------------------------------------
% Based on the generated sensitivity matrix "SA_matrix," correlation
% analysis can be conducted.
%--------------------------------------------------------------------------
% Rearrange both the sensitivity matrix and ESPM parameter sensitivity vector
% in descending order.
CA_use_matrix=SA_matrix(:,index_SA);
CA_use_SA=SA_ESPM(:,index_SA);

%--------------------------------------------------------------------------
% Calculate correlation matrix
% Details can be found in: https://ieeexplore.ieee.org/document/9195005 (Eq. 77 and Eq.78)
%--------------------------------------------------------------------------
for kk=1:1:13
    for jj=1:1:13
      Corr_matrix(kk,jj) =dot(CA_use_matrix(:,kk),CA_use_matrix(:,jj))/(CA_use_SA(kk)*CA_use_SA(jj)) ;
    end
end

Corr_matrix=fliplr(abs(Corr_matrix));

for kk=1:1:13
    Corr_matrix(13-kk+1,kk)=0;
end

CA_results{totkk}=Corr_matrix;


%--------------------------------------------------------------------------
% Based on the generated correlation matrix, conduct correlation analysis.
% For correlation analysis, we need to set a threshold.
% When the correlation value between two parameters exceeds this threshold,
% they are considered correlated parameters.
%--------------------------------------------------------------------------
% Define the vector for correlation analysis.
% In this example, we use 3 types of current data to run LSA and CA,
% therefore we need to set 3 threshold values.
threshold_CA=[0.95,0.95,0.99];
threshold=threshold_CA(totkk);




%% Design the parameter identification strategy based on 
%  sensitivity and correlation analysis results
function_iden_strategy(SA_flag,totkk,Corr_matrix,threshold,gsax_labels);
SA_flag=SA_flag+1;
end


%% Plot LSA and CA results
figure
fig_v1=[1,2,3];
fig_v2=[[4,7,10];[5,8,11];[6,9,12]];
title_v={'C/40 discharge','C/2 discharge','HPPC'};
for kk=1:1:3
%--------------------------------------------------------------------------
% Plot the sensitivity analysis results
%--------------------------------------------------------------------------
show_labels={'R_s_,_n','R_s_,_p','ε_n','ε_p','D_s_,_n','D_s_,_p','k_n','k_p','L_n','L_p','D_e','κ_e','R_s'};
gsax_labels=show_labels(SA_results{kk}(1,:));
x=1:1:13;
y=SA_results{kk}(2,:);

subplot(4,3,fig_v1(kk))
scatter(x, y, 'filled');
set(gca, 'YScale', 'log');
set(gca,'XTick',1:13,'XTickLabel',gsax_labels)
ylabel('Sensitivty');
title(title_v(kk))
grid on;
box on
hold on
set(gca,'linewidth',1,'fontsize',10,'fontname','Arial');


%--------------------------------------------------------------------------
% Plot the corrolation analysis results
%--------------------------------------------------------------------------
Corr_matrix=CA_results{kk};
threshold=threshold_CA(kk);
subplot(4,3,fig_v2(kk,:))
h = heatmap(Corr_matrix);
h.XDisplayLabels = fliplr(gsax_labels);
h.YDisplayLabels = gsax_labels;
cmap = colormap(h);
num_colors = size(cmap, 1);
threshold_color_index = round(threshold * num_colors);
% Set the color type for highly correlated parameters (exceeding the threshold).
new_cmap = cmap;
new_cmap(threshold_color_index:end, :) = repmat([1, 0, 0], num_colors - threshold_color_index + 1, 1);
colormap(h, new_cmap);
caxis([0 1]);
xlabel('ESPM parameters');
ylabel('ESPM parameters');
h.FontSize = 10; 
h.CellLabelColor = 'none'; 
set(gcf,'unit','centimeters','position',[13 5 20 20])
end

set(gcf,'unit','centimeters','position',[0 5 50 15])


% %--------------------------------------------------------------------------
% % Plot the identification stragety
% % 1. Model parameters identified in this step are shown in yellow.
% % 2. Model parameters obtained from the last step are shown in purple.
% % 3. Model parameters that are correlated with each other are shown in black.
% %--------------------------------------------------------------------------
% 
% fig = openfig('Iden_strategy1.fig');
% set(fig, 'Units', 'centimeters', 'Position', [0, 5, 15, 15]);
% title('C/40')
% 
% fig = openfig('Iden_strategy2.fig');
% set(fig, 'Units', 'centimeters', 'Position', [15, 5, 15, 15]);
% title('C/2')
% 
% fig = openfig('Iden_strategy3.fig');
% set(fig, 'Units', 'centimeters', 'Position', [30, 5, 15, 15]);
% title('HPPC')


%--------------------------------------------------------------------------
% Show identification stragety
% Profile 1 represents C/40 discharge
% Profile 2 represents C/2 discharge
% Profile 3 represents HPPC
%--------------------------------------------------------------------------
message = '==============================================================';
fprintf('%s\n', message);
message = 'This is the identification stragety based on LSA and CA.';
fprintf('%s\n', message);
message = '==============================================================';
fprintf('%s\n', message);
% Define multiple cell arrays
cellArray1=load('SA_re_1.mat');
cellArray1=cellArray1.iden_par.row;
cellArray2=load('SA_re_2.mat');
cellArray2=cellArray2.iden_par.row;
cellArray3=load('SA_re_3.mat');
cellArray3=cellArray3.iden_par.row;
% Combine them into a cell array of cell arrays
allCellArrays = {cellArray1, cellArray2, cellArray3};
% Loop through each cell array and display its content
for i = 1:length(allCellArrays)
    fprintf('Identified parameters using profile %d:\n', i);
    for j = 1:length(allCellArrays{i})
        fprintf('%s\n', allCellArrays{i}{j});
    end
end