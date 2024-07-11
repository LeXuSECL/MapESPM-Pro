% -------------------------------------------------------------------------
% ESPM_main.m
% -------------------------------------------------------------------------
% ESPM_main is the script for ESPM model simulation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MapESPM-Pro: Mass-Preserving Extended Single Particle Model with Parameter Optimization for Control-Oriented Applications

% Copyright (c) 2024 Mass-Preserving Extended Single Particle Model with Parameter Optimization for Control-Oriented Applications (MapESPM-Pro)
% MapESPM-Pro is freely distributed software under the MIT License 
% v1.0.0.: Released 07/2024 

% Written by Le Xu (lexu1209@stanford.edu)
% PI: Prof. Simona Onori (sonori@stanford.edu)
% Stanford Energy Control Group, Stanford University
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------------------------------------------------------------------
% The input for the ESPM function is:
% 1. ESPM parameter values: x
% 2. Time: t_data
% 3. Current: I_data
% 4. Initial SOC: SOC_IC
% 5. Ambient Temperature: T_amb
% 6. Data type: profile_flag

% The output of the ESPM function is:
% 1. Battery voltage: V_cell
% 2. Negative electrode SOC: soc_bulk_n
% 3. Positive electrode SOC: soc_bulk_p
% 4. Negative electrode concentration: cs_n
% 5. Positive electrode concentration: cs_p
% 6. Open-circuit voltage: V_oc
% 7. Parameters used for simulation: param
% 8. Positive electrode potential: ocp_p
% 9. Negative electrode potential: ocp_n
% 10. Positive electrode overpotential: eta_p
% 11. Negative electrode overpotential: eta_n
% 12. Electrolyte concentration: ce
% 13. Electrolyte overpotential: eta_ele
% -------------------------------------------------------------------------

function [V_cell, R_l, T_core, T_surf,soc_bulk_n, soc_bulk_p, cs_n, cs_p,...
          V_oc,param,ocp_p,ocp_n,eta_p,eta_n,ce,eta_ele] = ESPM_main(x, t_data, I_data, ...
   SOC_IC,T_amb,profile_flag) 

% -------------------------------------------------------------------------
% SPECIFY: Finite Volume Method (FVM) Discretization Settings
% -------------------------------------------------------------------------
param.Nr = 21;      % Number of radial discretization grids in ESPM
param.Nx_n = 10;    % Number of cartesian discretization grids in ESPM
param.Nx_s = 10;    % Number of cartesian discretization grids in ESPM
param.Nx_p = 10;    % Number of cartesian discretization grids in ESPM


% -------------------------------------------------------------------------
% Load input for ESPM
% -------------------------------------------------------------------------
param.profile_flag=profile_flag;
param.t_data = t_data;
param.I_data = I_data;                        

% -------------------------------------------------------------------------
% Initialize all model parameters
% -------------------------------------------------------------------------
run ModelParameters.m

% -------------------------------------------------------------------------
% SPECIFY: initial conditions of all cells
% -------------------------------------------------------------------------
% Set Ambient Temperature
param.T_amb = 273 + T_amb; %[K]

% Set Initial SOC
SOC_cells = SOC_IC;  

% Set Initial lithium concentration in solid-phase
[cs_initial,csn0,csp0] = conc_initial_sd(SOC_cells,param);

% Set Initial lithium concentration in electrolyte
param.ce_states = param.Nx_n + param.Nx_s + param.Nx_p;
ce_initial = param.ce0*ones(param.ce_states,1);

% Initial temperature for all cells
param.T_states = 2; % Factor of 2 because of two-state model


% Combine all initial conditions to x_initial
x_initial = [cs_initial;ce_initial];     

% -------------------------------------------------------------------------
% Run ESPM Simulation
% -------------------------------------------------------------------------
 [V_cell, R_l, T_core, T_surf,soc_bulk_n, soc_bulk_p, cs_n, cs_p,...
          V_oc,param,ocp_p,ocp_n,eta_p,eta_n,ce,eta_ele] = Module_Run_casadi(x_initial,param);
end



function [cs_initial, csn0, csp0] = conc_initial_sd(SOC_cells, param)
% According to definiation, SOC in electrochemical model is calculated as
%        θp0-θp,bulk
%SOCp= ----------------   
%         θp0-θp100 

%        θn,bulk-θn,0
%SOCn= ----------------
%         θn100-θn0 

csn0 = (SOC_cells*(param.theta100_n - param.theta0_n) + param.theta0_n).*param.c_n_max ;
csp0 = (param.theta0_p - SOC_cells*(param.theta0_p - param.theta100_p)).*param.c_p_max; 

% Each electrode is divided into N control volumes with same inital values
csn_initial = csn0*ones((param.Nr-1),1);
csp_initial = csp0*ones((param.Nr-1),1);

cs_initial = [csn_initial;csp_initial];
end


function [V_cell, R_l, T_core, T_surf,soc_bulk_n, soc_bulk_p, cs_n, cs_p,...
          V_oc,param,ocp_p,ocp_n,eta_p,eta_n,ce,eta_ele] = Module_Run_casadi(x_initial,param)
% Load the solver      
import casadi.*      
global int_method model_type


param.TIMEOLD = datetime('now','Format','HHmmss');

%% Define all the governing equations for ESPM model
totN=size(x_initial,1);
% Define all the ESPM variables
% 1 csn  2 csp  3 ce
sundials_u = SX.sym('sundials_u',totN);
input_crt=SX.sym('input_crt');
fun_ode=Module_ode_casadi(sundials_u,input_crt, param);
dae = struct('x',sundials_u,'p',input_crt,'ode',[fun_ode],'alg',[]);
%%

% Load the initial values
u0ini=x_initial;
u0ini_casadi=u0ini;
slice_index=find(abs(param.I_data(2:end)-param.I_data(1:end-1))>1e-3);
slice_index=[0,slice_index,size(param.I_data,2)];
tot_slice=size(slice_index,2);
rextot=[];

%% Casadi solver settings
opts_defult= struct('tf',1,'abstol',1e-3,'reltol',1e-3,'linear_solver','csparse',...
                'newton_scheme','gmres');  
F1_defult = integrator('F', 'cvodes', dae,opts_defult);

% Solve the model using Casadi
for jj=1:1:tot_slice-1
    Iapp=param.I_data(slice_index(jj)+1:slice_index(jj+1));
    rex=[];
    if size(Iapp,2)<10   % In this mode, Casadi will solve step by step
                for kk=1:1:size(Iapp,2)
                iapp=Iapp(kk);
                sol = F1_defult('x0',u0ini_casadi,'p',iapp);
                u0ini_casadi=full(sol.xf);
                rex(:,kk)=full(sol.xf);
                end
                 
    else            % In this mode, Casadi will solve all at once        
                tgrid=1:1:size(Iapp,2)+1;
                opts = struct('grid',tgrid,'output_t0',false,'abstol',1.0e-3,'reltol',1e-3,'linear_solver','csparse',...
                'newton_scheme','gmres','show_eval_warnings',false); 
                F1 = integrator('F', 'cvodes', dae,opts);
                iapp=mean(Iapp);
                sol = F1('x0',u0ini_casadi,'p',iapp);              
                rex=full(sol.xf);
                u0ini_casadi=rex(:,end);
    end
    rextot=[rextot,rex];        
end
rextot=[u0ini,rextot(:,1:end-1)]; % This is because casadi didn't save the initial values
x_out=rextot;


%% calculate surface concentration 
% Using FVM, we need to further calculate surface concentraion 
% Because in FVM, we only obtain average concentraion at each CV
% More details can be found in here:
% https://iopscience.iop.org/article/10.1149/1945-7111/ad1293

% =============== Average concentraion at each CV in FVM method
cs = x_out(1:(param.Nr-1)*2,:);            % All solid concentrations
cs_n = cs(1:(param.Nr-1),:);               % Negative electrode Concentrations
cs_p = cs((param.Nr-1)+1:end,:);           % Positive electrode Concentrations

% ========= To calculate surface concentraion, we need interpolation method
% Here we provide two interpolation methods
% int_method ==1 The linear interpolation method
% int_method ==2 The Hermite interpolation method (We proposed this)
% More details for int_method==2 can be found in here:
% https://iopscience.iop.org/article/10.1149/1945-7111/ad1293
if int_method==2 
% First, calculate the surface concentraion for positive electrode
dr_p=param.Rs_p/(param.Nr-1);
rawx=linspace(dr_p/2,param.Rs_p-dr_p/2,param.Nr-1);
rawy=cs_p;
for k=1:1:size(rawy,2)
addy(k)=mypchip(rawx(end-2:end),rawy(end-2:end,k),param.Rs_p); % mypchip is the Hermite interpolation function
end
new_cssurp2=addy;

% Then, calculate the surface concentraion for negative electrode
dr_n=param.Rs_n/(param.Nr-1);
rawx=linspace(dr_n/2,param.Rs_n-dr_n/2,param.Nr-1);
rawy=cs_n;
for k=1:1:size(rawy,2)
addy(k)=mypchip(rawx(end-2:end),rawy(end-2:end,k),param.Rs_n); % mypchip is the Hermite interpolation function
end
new_cssurn2=addy;
elseif int_method==1
new_cssurp2=(3*cs_p(end,:)-cs_p(end-1,:))/2;
new_cssurn2=(3*cs_n(end,:)-cs_n(end-1,:))/2;    
end


% ===== Obtain electrolyte concentraion
ce=x_out(end-param.ce_states+1:end,:);

% Load cell temperature
T_cell=param.T_amb*ones(2,1)*ones(1,size(x_out,2));
T_core(1,:) = T_cell(1,:);
T_surf(1,:) = T_cell(2,:);

% Calculation θn and θp according to:
%       csurf,n
%θn= ----------
%       csmax,n       

%       csurf,p
%θp= ----------
%       csmax,p      

theta_surf_n(1,:) = new_cssurn2(:)/param.c_n_max;
theta_surf_p(1,:) = new_cssurp2(:)/param.c_p_max;   



% Load kn and kp
kn = param.kn_ref;
kp = param.kp_ref;


%% Becasue we have hystersis, we need to use different Up and Un 
% Find the timesteps when battery is been charged or discharged
index_cc=find(param.I_data<0);
index_dc=find(param.I_data>0);
index_rest=find(param.I_data==0);

ocp_p(index_cc)=U_p_cc(theta_surf_p(index_cc));
ocp_p(index_dc)=U_p_dis(theta_surf_p(index_dc));
ocp_p(index_rest)=(U_p_dis(theta_surf_p(index_rest))+U_p_dis(theta_surf_p(index_rest)))/2;

ocp_n(index_cc)=U_n_cc(theta_surf_n(index_cc));
ocp_n(index_dc)=U_n_dis(theta_surf_n(index_dc));
ocp_n(index_rest)=(U_n_dis(theta_surf_n(index_rest))+U_n_dis(theta_surf_n(index_rest)))/2;
%%

% ======Calculate overpotential
% 1 Overpotential from positive electrode
eta_p(1,:) = eta_cathode(theta_surf_p(1,:), ce(1:param.Nx_p,:), T_core(1,:), param.I_data(1,:), param, kp);
% 2 Overpotential from negative electrode
eta_n(1,:) = eta_anode(theta_surf_n(1,:), ce(end-param.Nx_n+1:end,:), T_core(1,:), param.I_data(1,:), param, kn);  
% 3 Overpotential from electrolyte
eta_ele(1,:)=eta_electrolyte(ce(:,:),param.I_data(1,:),param);

% Calculate open-circuit voltage
% OCV=Up - Un
V_oc = ocp_p - ocp_n;

%% Calculate SOC
% 1 Calculate the volume average SOC for each electrode
csp_ave=volume_ave_con_fvm(param,cs_p,0);
csn_ave=volume_ave_con_fvm(param,cs_n,1);
% 2 Calculate bulk stochiomitery value for each electrode
xbulkn=csn_ave/param.c_n_max;
xbulkp=csp_ave/param.c_p_max;
% 3 Calculate bulk SOC for each electrode
soc_bulk_p=(param.theta0_p-xbulkp)/(param.theta0_p - param.theta100_p);
soc_bulk_n=(xbulkn - param.theta0_n)/(param.theta100_n - param.theta0_n);

% Calculate terminal voltage for the battery
[V_cell(1,:), R_l(1,:)] = V_calculation(ocp_p(1,:),ocp_n(1,:),eta_p(1,:),eta_n(1,:),eta_ele(1,:),...
    param.I_data(1,:),param,model_type);
end


% -------------------------------------------------------------------------
% Obtain ODEs for solid-phase and electrolyte
% -------------------------------------------------------------------------
function [dx_dt] = Module_ode_casadi(x_in,input_crt, param)
import casadi.*
% Discritize solid phase diffusion equation based on FVM
dcs_dt=ODE_solid_phase_FVM(x_in(1:end-param.ce_states),input_crt, param);
% Discritize electrolyte diffusion equation based on FVM
dce_dt=ODE_electrolyte_FVM(x_in(end-param.ce_states+1:end),input_crt, param);

dx_dt = [dcs_dt;dce_dt];
end

% -------------------------------------------------------------------------
% Function to Calculate Negative Electrode Overpotential
% -------------------------------------------------------------------------
function eta_n = eta_anode(x, ce_n, T_core, input_crt, param, kn)
ce_n_avg=mean(ce_n);
i0_n= kn*param.F.*(ce_n_avg.^0.5).*...
    ((param.c_n_max.*x).^0.5).*...
    ((param.c_n_max-param.c_n_max.*x).^0.5); 
eta_n = asinh(input_crt./(2*param.A*(param.a_sn)*...
    param.Ln.*i0_n)).*(param.Rg.*T_core)./...
    (param.F*param.alpha_cell);
end

% -------------------------------------------------------------------------
% Function to Calculate Positive Electrode Overpotential
% -------------------------------------------------------------------------
function eta_p = eta_cathode(x, ce_p, T_core, input_crt, param, kp)
ce_p_avg=mean(ce_p);
i0_p= kp*param.F.*(ce_p_avg.^0.5).*...
    ((param.c_p_max.*x).^0.5).*...
    ((param.c_p_max-param.c_p_max.*x).^0.5);   
eta_p = asinh(-input_crt./(2*param.A*(param.a_sp)*...
    param.Lp.*i0_p)).*(param.Rg.*T_core)./...
    (param.F*param.alpha_cell);
end

% -------------------------------------------------------------------------
% Function to calculate electrolyte overpotential
% -------------------------------------------------------------------------
function etae = eta_electrolyte(ce,Iinput,param)
eps_el_p=param.eps_el_p;
eps_el_s=param.eps_el_s;
eps_el_n=param.eps_el_n_ref;

kep=param.kep;
kes=param.kes;
ken=param.ken;

keffp=kep*eps_el_p^param.brugg_p;
keffs=kes*eps_el_s^param.brugg_s;
keffn=ken*eps_el_n^param.brugg_n;

Lp=param.Lp;
Ls=param.Ls;
Ln=param.Ln;

etae1=2*param.Rg*param.T_amb*(1-param.t0)/param.F*log(ce(1,:)./ce(end,:));
etae2=Iinput/(2*param.A)*(Lp/keffp+2*Ls/keffs+Ln/keffn);
etae=etae1-etae2;
end

% -------------------------------------------------------------------------
% Function for Un charge
% -------------------------------------------------------------------------
function y = U_n_cc(x)
% Un charge as a function of SOC
% Below are the parameters of the Un charge curve
par=[1056906.89598384,0,-2689073.86515939,0,0,-4621784.83877322,509132.488750815,5611545.23157032,3032506.44367608,-2110204.58197130,-3959976.54182831,-160624.788244324,1854791.95269040,266589.673231300,-345689.450777555,-42426.6577404344,15322.0182869209,1.04719443819735];
a0=par(1);
a1=par(2);
b1=par(3);
a2=par(4);
b2=par(5);
a3=par(6);
b3=par(7);
a4=par(8);
b4=par(9);
a5=par(10);
b5=par(11);
a6=par(12);
b6=par(13);
a7=par(14);
b7=par(15);
a8=par(16);
b8=par(17);
w=par(18);
% x is the SOC, given SOC as input, the output is OCP
y=a0 + a1*cos(x*w) + b1*sin(x*w) + ...
a2*cos(2*x*w) + b2*sin(2*x*w) + a3*cos(3*x*w) + b3*sin(3*x*w) + ...
a4*cos(4*x*w) + b4*sin(4*x*w) + a5*cos(5*x*w) + b5*sin(5*x*w) + ...
a6*cos(6*x*w) + b6*sin(6*x*w) + a7*cos(7*x*w) + b7*sin(7*x*w) + ...
a8*cos(8*x*w) + b8*sin(8*x*w);
end

% -------------------------------------------------------------------------
% Function for Un discharge
% -------------------------------------------------------------------------
function y = U_n_dis(x)
% Un discharge as a function of SOC
% Below are the parameters of the Un discharge curve
par=[-77379.6258610637,0,211536.280133149,0,0,394664.275692445,-78146.7574323678,-546030.847477876,-218492.590911161,261424.418592925,361318.557117052,-13252.4306514239,-200845.609381728,-24542.0560032409,45489.2190270259,5117.79196290630,-2887.95721989976,1.04719578297250];
a0=par(1);
a1=par(2);
b1=par(3);
a2=par(4);
b2=par(5);
a3=par(6);
b3=par(7);
a4=par(8);
b4=par(9);
a5=par(10);
b5=par(11);
a6=par(12);
b6=par(13);
a7=par(14);
b7=par(15);
a8=par(16);
b8=par(17);
w=par(18);
% x is the SOC, given SOC as input, the output is OCP
y=a0 + a1*cos(x*w) + b1*sin(x*w) + ...
a2*cos(2*x*w) + b2*sin(2*x*w) + a3*cos(3*x*w) + b3*sin(3*x*w) + ...
a4*cos(4*x*w) + b4*sin(4*x*w) + a5*cos(5*x*w) + b5*sin(5*x*w) + ...
a6*cos(6*x*w) + b6*sin(6*x*w) + a7*cos(7*x*w) + b7*sin(7*x*w) + ...
a8*cos(8*x*w) + b8*sin(8*x*w);
end

% -------------------------------------------------------------------------
% Function for Up charge
% -------------------------------------------------------------------------
function y = U_p_cc(x)
% Up charge as a function of SOC
% Below are the parameters of the Up charge curve
par=[36439928.6635994,-41721043.6266609,-25431018.9200839,0,0,-2326719.38249511,32445226.6017429,19643144.6396902,-27559278.1451687,-17637319.7219740,7534682.97551583,6599125.42991595,959131.376944491,-1036700.98011868,-860374.830484993,39589.4839541989,115858.079543840,1.04719756661260];
a0=par(1);
a1=par(2);
b1=par(3);
a2=par(4);
b2=par(5);
a3=par(6);
b3=par(7);
a4=par(8);
b4=par(9);
a5=par(10);
b5=par(11);
a6=par(12);
b6=par(13);
a7=par(14);
b7=par(15);
a8=par(16);
b8=par(17);
w=par(18);
% x is the SOC, given SOC as input, the output is OCP
y=a0 + a1*cos(x*w) + b1*sin(x*w) + ...
a2*cos(2*x*w) + b2*sin(2*x*w) + a3*cos(3*x*w) + b3*sin(3*x*w) + ...
a4*cos(4*x*w) + b4*sin(4*x*w) + a5*cos(5*x*w) + b5*sin(5*x*w) + ...
a6*cos(6*x*w) + b6*sin(6*x*w) + a7*cos(7*x*w) + b7*sin(7*x*w) + ...
a8*cos(8*x*w) + b8*sin(8*x*w);
end

% -------------------------------------------------------------------------
% Function for Up discharge
% -------------------------------------------------------------------------
function y = U_p_dis(x)
% Up discharge as a function of SOC
% Below are the parameters of the Up discharge curve
par=[-39033981.1321803,47127267.0820635,22798594.8315349,0,0,-7585030.20786143,-34088437.8332242,-8361427.70351725,35422374.5930115,13040230.4868199,-16038555.1349497,-6528177.95774521,3051057.71037441,1460174.08797285,20197.2134023230,-119050.203844771,-59717.2026753066,1.04719756256240];
a0=par(1);
a1=par(2);
b1=par(3);
a2=par(4);
b2=par(5);
a3=par(6);
b3=par(7);
a4=par(8);
b4=par(9);
a5=par(10);
b5=par(11);
a6=par(12);
b6=par(13);
a7=par(14);
b7=par(15);
a8=par(16);
b8=par(17);
w=par(18);
% x is the SOC, given SOC as input, the output is OCP
y=a0 + a1*cos(x*w) + b1*sin(x*w) + ...
a2*cos(2*x*w) + b2*sin(2*x*w) + a3*cos(3*x*w) + b3*sin(3*x*w) + ...
a4*cos(4*x*w) + b4*sin(4*x*w) + a5*cos(5*x*w) + b5*sin(5*x*w) + ...
a6*cos(6*x*w) + b6*sin(6*x*w) + a7*cos(7*x*w) + b7*sin(7*x*w) + ...
a8*cos(8*x*w) + b8*sin(8*x*w);
end

% -------------------------------------------------------------------------
% This is the Hermite Interpolation Method
% -------------------------------------------------------------------------
function yy = mypchip(x,y,xx)
% More details can be found in here:
% https://iopscience.iop.org/article/10.1149/1945-7111/ad1293
y1=y(2);
y2=y(3);
x1=x(2);
x2=x(3);
dy1=(y(3)-y(1))/(x(3)-x(1));
dy2=(3*y(3)-4*y(2)+y(1))/(x(3)-x(1));
term1=(1+2*(xx-x1)/(x2-x1))*((xx-x2)/(x1-x2))^2*y1;
term2=(1+2*(xx-x2)/(x1-x2))*((xx-x1)/(x2-x1))^2*y2;
term3=(xx-x1)*((xx-x2)/(x1-x2))^2*dy1;
term4=(xx-x2)*((xx-x1)/(x2-x1))^2*dy2;
yy=term1+term2+term3+term4;
end

% -------------------------------------------------------------------------
% Electrolyte Governing Equation Discretized by FVM
% -------------------------------------------------------------------------
function dy = ODE_electrolyte_FVM(u, input_crt,data)
import casadi.*

% Load parameter settings
L1=data.Lp;
L2=data.Ls;
L3=data.Ln;
N1=data.Nx_p;
N2=data.Nx_s;
N3=data.Nx_n;
dx1=L1/N1;
dx2=L2/N2;
dx3=L3/N3;
totN=N1+N2+N3;
Iinput=input_crt;
A=data.A;
bruggp=data.brugg_p;
bruggs=data.brugg_s;
bruggn=data.brugg_n;
eps_el_p=data.eps_el_p;
eps_el_s=data.eps_el_s;
eps_el_n=data.eps_el_n_ref;
lp=L1;
ln=L3;
ap=data.a_sp;
an=data.a_sn;
tp=data.t0;
F=data.F;
j_p = -Iinput/(ap*F*lp*A);
j_n = Iinput/(an*F*ln*A);


De=data.De;
De_p=De*ones(N1,1);
De_s=De*ones(N2,1);
De_n=De*ones(N3,1);
Deff_p=De_p*eps_el_p^bruggp; 
Deff_s=De_s*eps_el_s^bruggs; 
Deff_n=De_n*eps_el_n^bruggn; 

Deff_p_m05=Deff_p(1:end-1);
Deff_p_p05=Deff_p(2:end);

Deff_s_m05=Deff_s(1:end-1);
Deff_s_p05=Deff_s(2:end);

Deff_n_m05=Deff_n(1:end-1);
Deff_n_p05=Deff_n(2:end);

% -------------------------------------------------------------------------
% Using harmonic mean (HM) 
% More details can be found here: (LIONSIMBA Paper) Page7
% https://iopscience.iop.org/article/10.1149/2.0291607jes 
% -------------------------------------------------------------------------
Deff_p_15=2*Deff_p_m05.*Deff_p_p05./(Deff_p_m05+Deff_p_p05);
Deff_s_15=2*Deff_s_m05.*Deff_s_p05./(Deff_s_m05+Deff_s_p05);
Deff_n_15=2*Deff_n_m05.*Deff_n_p05./(Deff_n_m05+Deff_n_p05);

% -------------------------------------------------------------------------
% Harmonic mean for the Positive/Separator and Negtaive/Separator interface
% More details can be found here: (LIONSIMBA Paper) Page7
% https://iopscience.iop.org/article/10.1149/2.0291607jes 
% -------------------------------------------------------------------------
beta_ps=dx1/(dx1+dx2);
beta_sn=dx2/(dx2+dx3);
% Define Interface
% Positive and Separator interface
Deff_ps=Deff_p(end)*Deff_s(1)/(beta_ps*Deff_s(1)+(1-beta_ps)*Deff_p(end));
% Negative and Separator interface
Deff_sn=Deff_s(end)*Deff_n(1)/(beta_sn*Deff_n(1)+(1-beta_sn)*Deff_s(end));


% -------------------------------------------------------------------------
% Define ODE for Positive, Separator and Negative
% -------------------------------------------------------------------------
y=u;
dy=SX.sym('dy',totN);

% Positive Region
dy(1)=(Deff_p_15(1)*(y(2)-y(1))/(dx1)^2+ap*(1-tp)*j_p)/eps_el_p;
for kk=2:1:N1-1
 dy(kk)=(Deff_p_15(kk)*(y(kk+1)-y(kk))/(dx1)^2-Deff_p_15(kk-1)*(y(kk)-y(kk-1))/(dx1)^2+ap*(1-tp)*j_p)/eps_el_p;       
end
dy(N1)=(Deff_ps*(y(N1+1)-y(N1))/((dx1+dx2)*(dx1)/2)-Deff_p_15(N1-1)*(y(N1)-y(N1-1))/(dx1)^2+ap*(1-tp)*j_p)/eps_el_p;
% Separator Region
dy(N1+1)=(Deff_s_15(1)*(y(N1+2)-y(N1+1))/(dx2)^2-Deff_ps*(y(N1+1)-y(N1))/((dx1+dx2)*(dx2)/2))/eps_el_s;
for kk=N1+2:1:N1+N2-1
dy(kk)=(Deff_s_15(kk-N1)*(y(kk+1)-y(kk))/(dx2)^2-Deff_s_15(kk-N1-1)*(y(kk)-y(kk-1))/(dx2)^2)/eps_el_s;   
end
dy(N1+N2)=(Deff_sn*(y(N1+N2+1)-y(N1+N2))/((dx2+dx3)*(dx2)/2)-Deff_s_15(N2-1)*(y(N1+N2)-y(N1+N2-1))/(dx2)^2)/eps_el_s;
% Negative Region
dy(N1+N2+1)=(Deff_n_15(1)*(y(N1+N2+2)-y(N1+N2+1))/(dx3)^2-Deff_sn*(y(N1+N2+1)-y(N1+N2))/((dx2+dx3)*(dx3)/2)+an*(1-tp)*j_n)/eps_el_n;
for kk=N1+N2+2:1:N1+N2+N3-1
dy(kk)=(Deff_n_15(kk-N1-N2)*(y(kk+1)-y(kk))/(dx3)^2-Deff_n_15(kk-N1-N2-1)*(y(kk)-y(kk-1))/(dx3)^2+an*(1-tp)*j_n)/eps_el_n;   
end
dy(N1+N2+N3)=(-Deff_n_15(N3-1)*(y(N1+N2+N3)-y(N1+N2+N3-1))/(dx3)^2+an*(1-tp)*j_n)/eps_el_n; 
end


% -------------------------------------------------------------------------
% Solid-phase Diffusion Governing Equation Discretized by FVM
% -------------------------------------------------------------------------
function dudt = ODE_solid_phase_FVM(u_in, input_crt,param)
% Load parameter settings
Ncv=param.Nr-1;
u_in=u_in(1:2*param.Nr-2);
us_n=u_in(1:param.Nr-1);
us_p=u_in(param.Nr:end);
dr_n=param.Rs_n/(Ncv);
dr_p=param.Rs_p/(Ncv);
ap=param.a_sp;
an=param.a_sn;
F=param.F;
Lp=param.Lp;
Ln=param.Ln;

% -------------------------------------------------------------------------
% Define ODE for Positive and Negative electrodes
% More details can be found in here:
% https://iopscience.iop.org/article/10.1149/1945-7111/ad1293
% -------------------------------------------------------------------------
% Positive Electrode
n3p=-3*(Ncv)^2*(-1)/(dr_p*(3*(Ncv)^2-3*(Ncv)+1)*(F*ap*Lp*param.A));
Ap=[];
for kk=1:1:Ncv   
if kk==1
    Ap(kk,1)=-1*((3*kk^2)/(4*(3*kk^2 - 3*kk + 1)));
    Ap(kk,2)= 1*((3*kk^2)/(4*(3*kk^2 - 3*kk + 1)));
elseif kk>1&&kk<Ncv
    Ap(kk,kk-1)=1*((3*(kk - 1)^2)/(4*(3*kk^2 - 3*kk + 1)));
    Ap(kk,kk)= -1*((6*kk^2 - 6*kk + 3)/(4*(3*kk^2 - 3*kk + 1)));
    Ap(kk,kk+1)= 1*((3*kk^2)/(4*(3*kk^2 - 3*kk + 1)));
elseif kk==Ncv
    Ap(kk,kk-1)= 1*((3*((Ncv) - 1)^2)/(4*(3*(Ncv)^2 - 3*(Ncv) + 1)));
    Ap(kk,kk)= -1*((3*((Ncv) - 1)^2)/(4*(3*(Ncv)^2 - 3*(Ncv) + 1)));
end
end

Bp=zeros(size(us_p,1),1);
for kk=1:1:Ncv
if kk~=Ncv
    Bp(kk)=0;
else
    Bp(kk)=n3p;
end
end

% Negative Electrode
n3n=-3*(Ncv)^2*(1)/(dr_n*(3*(Ncv)^2-3*(Ncv)+1)*(F*an*Ln*param.A));
An=[];
for kk=1:1:Ncv    
if kk==1
    An(kk,1)=-1*((3*kk^2)/(4*(3*kk^2 - 3*kk + 1)));
    An(kk,2)= 1*((3*kk^2)/(4*(3*kk^2 - 3*kk + 1)));
elseif kk>1&&kk<Ncv
    An(kk,kk-1)=1*((3*(kk - 1)^2)/(4*(3*kk^2 - 3*kk + 1)));
    An(kk,kk)= -1*((6*kk^2 - 6*kk + 3)/(4*(3*kk^2 - 3*kk + 1)));
    An(kk,kk+1)= 1*((3*kk^2)/(4*(3*kk^2 - 3*kk + 1)));
elseif kk==Ncv
    An(kk,kk-1)= 1*((3*((Ncv) - 1)^2)/(4*(3*(Ncv)^2 - 3*(Ncv) + 1)));
    An(kk,kk)= -1*((3*((Ncv) - 1)^2)/(4*(3*(Ncv)^2 - 3*(Ncv) + 1)));
end
end

Bn=zeros(size(us_n,1),1);
for kk=1:1:Ncv
if kk~=Ncv
    Bn(kk)=0;
else
    Bn(kk)=n3n;
end
end

dudt_un=An*4*param.Dsn_ref/dr_n^2*us_n+Bn*input_crt;
dudt_up=Ap*4*param.Dsp_ref/dr_p^2*us_p+Bp*input_crt;
dudt=[dudt_un;dudt_up];
end

% -------------------------------------------------------------------------
% Function to calculate battery voltage
% -------------------------------------------------------------------------
function [V_cell, R_l] = V_calculation(ocp_p,ocp_n,eta_p,eta_n,eta_e,I, param, model_type)
R_l=param.R_l;
if model_type==1  % SPM model
V_cell = ocp_p - ocp_n + eta_p - eta_n-I.*R_l;
elseif model_type==2  % ESPM model
V_cell = ocp_p - ocp_n + eta_p - eta_n+eta_e-I.*R_l;
end
end

% -------------------------------------------------------------------------
% Function to calculate volume-average solid-phase concentraion
% -------------------------------------------------------------------------
function csi_ave = volume_ave_con_fvm(param,cs_i,flag)
% The flag indicate positive electrode (0) or negative electrode (1)
if flag==0   % Positive electrode
cs_p=cs_i;    
dr_p=param.Rs_p/(param.Nr-1);
r_p=linspace(dr_p,param.Rs_p,param.Nr-1);
r_p=[0,r_p];
for kk=2:1:length(r_p)
avecspi(kk-1,:)=cs_p(kk-1,:)*4/3*pi*(r_p(kk)^3-r_p(kk-1)^3);
end
csi_ave=sum(avecspi)/(4/3*pi*(param.Rs_p^3));

elseif flag==1  % Negative electrode
cs_n=cs_i;   
dr_n=param.Rs_n/(param.Nr-1);
r_n=linspace(dr_n,param.Rs_n,param.Nr-1);
r_n=[0,r_n];
for kk=2:1:length(r_n)
avecsni(kk-1,:)=cs_n(kk-1,:)*4/3*pi*(r_n(kk)^3-r_n(kk-1)^3);
end
csi_ave=sum(avecsni)/(4/3*pi*(param.Rs_n^3));       
end

end

