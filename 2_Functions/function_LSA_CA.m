function [x,y,gsax_labels,SA_matrix,index_SA,SA_ESPM] = function_LSA_CA(SA_delta,data,profile_flag)

for SAkk=1:1:14
SAkk
espm_par=load('check.mat');
espm_par=espm_par.x_opt;
x_opt=espm_par;

if SAkk==1
   x_opt=x_opt; 
elseif SAkk==2
    x_opt(1)=-log10(10^(-x_opt(1))*(1-SA_delta));
elseif SAkk==3
    x_opt(2)=-log10(10^(-x_opt(2))*(1-SA_delta));  
elseif SAkk==4
    x_opt(3)=x_opt(3)*(1-SA_delta);
elseif SAkk==5
    x_opt(4)=x_opt(4)*(1-SA_delta);
elseif SAkk==6
    x_opt(9)=-log10(10^(-x_opt(9))*(1-SA_delta));
elseif SAkk==7
    x_opt(10)=-log10(10^(-x_opt(10))*(1-SA_delta));
elseif SAkk==8
    x_opt(11) = x_opt(11)-log10(1-SA_delta);
elseif SAkk==9
    x_opt(12) = x_opt(12)-log10(1-SA_delta);
elseif SAkk==10
    x_opt(17)=log10(10^(x_opt(17))*(1-SA_delta));
elseif SAkk==11
    x_opt(18)=log10(10^(x_opt(18))*(1-SA_delta));
elseif SAkk==12
    x_opt(20)=log10(10^(x_opt(20))*(1-SA_delta));
elseif SAkk==13
    x_opt(22)=x_opt(22)*(1-SA_delta);
elseif SAkk==14
    x_opt(28)=log10(10^(x_opt(28))*(1-SA_delta));
end
I_data=data.I_data;
t_data=data.t_data;
SOC_IC=data.SOC_IC;
T_amb=23;

int_method=data.int_method;
model_type=data.model_type;
battery_type=data.battery_type;

[V_cell, R_l, T_core, T_surf,soc_bulk_n, soc_bulk_p, cs_n, cs_p,...
          V_oc,param,ocp_p,ocp_n,eta_p,eta_n,ce,eta_ele] = ESPM_main(x_opt,t_data,I_data,SOC_IC,T_amb,profile_flag,int_method,model_type,battery_type);
      
SA_tot.V{SAkk}=V_cell';
SA_tot.socp{SAkk}=soc_bulk_p';
SA_tot.socn{SAkk}=soc_bulk_n';

SA_par(1)=10^(-x_opt(1));
SA_par(2)=10^(-x_opt(2));
SA_par(3)=x_opt(3);
SA_par(4)=x_opt(4);
SA_par(5)=10^(-x_opt(9));
SA_par(6)=10^(-x_opt(10));
SA_par(7)=10^(-x_opt(11))/param.F; 
SA_par(8)=10^(-x_opt(12))/param.F; 
SA_par(9)=10^(x_opt(17));
SA_par(10)=10^(x_opt(18));
SA_par(11)=10^(x_opt(20));
SA_par(12)=x_opt(22);
SA_par(13)=10^(x_opt(28));
SA_tot.par(:,SAkk)=SA_par;
end

%% After generating the values for the sensitivity matrix,
% the second step is to run sensitivity analysis.

dvdp=[];
nordvdq=[];
dsocpdp=[];
nordsocpdp=[];
dsocndp=[];
nordsocndp=[];
for kk=2:1:14
dvdp(:,kk-1)=(SA_tot.V{kk}- SA_tot.V{1})/(SA_tot.par(kk-1,kk)- SA_tot.par(kk-1,1));
nordvdq(:,kk-1)=dvdp(:,kk-1)*SA_tot.par(kk-1,1)/mean(SA_tot.V{1});
dsocpdp(:,kk-1)=(SA_tot.socp{kk}- SA_tot.socp{1})/(SA_tot.par(kk-1,kk)- SA_tot.par(kk-1,1));
nordsocpdp(:,kk-1)=dsocpdp(:,kk-1)*SA_tot.par(kk-1,1)/mean(SA_tot.V{1});
dsocndp(:,kk-1)=(SA_tot.socn{kk}- SA_tot.socn{1})/(SA_tot.par(kk-1,kk)- SA_tot.par(kk-1,1));
nordsocndp(:,kk-1)=dsocndp(:,kk-1)*SA_tot.par(kk-1,1)/mean(SA_tot.V{1});
end

% In sensitivity analysis, both voltage and 
% the State of Charge of the positive (SOCp) and negative (SOCn) electrodes are used.
SA_matrix=nordvdq+nordsocpdp+nordsocndp;
SA_ESPM=vecnorm(SA_matrix, 2, 1); % The Euclidean norm of each column represents the sensitivity of the corresponding parameter.

gsax_labels={'R_s_,_n','R_s_,_p','ε_n','ε_p','D_s_,_n','D_s_,_p','k_n','k_p','L_n','L_p','D_e','κ_e','R_s'};
x = linspace(1, 13, 13);
y = SA_ESPM;
[~,index_SA]=sort(y,'descend');
y=y(index_SA);
gsax_labels=gsax_labels(index_SA);
end

