function [x_opt, fval, history] = pso_for_step2(swarm_size, max_stalls, Initial_position, Lower_bound, Upper_bound,data)
model_type=data.model_type;
int_method=data.int_method;
battery_type=data.battery_type;
x_opt_raw=data.x_opt;
I_data=data.I_data;
t_data=data.t_data;
SOC_IC=data.SOC_IC;

refV=data.refV;
refsocp=data.refsocp;
refsocn=data.refsocn;
profile_flag=data.profile_flag;


obj_err=data.obj;

    history = [];

    
    options = optimoptions('particleswarm',...
                       'SwarmSize',swarm_size,...
                       'MaxStallIterations',max_stalls,...
                       'PlotFcn','pswplotbestf',...
                       'FunctionTolerance',1e-4,...
                       'ObjectiveLimit',obj_err,...
                       'InitialSwarmMatrix',Initial_position,...
                       'Display','iter',...
                       'SocialAdjustmentWeight',2,...
                       'SelfAdjustmentWeight',0.2, ...
                       'UseParallel',true, ...
                       'OutputFcn', @myoutput);

    [x_opt,fval,exitflag,output] = particleswarm(@(x) ESPM_fit_geom1(x),...
        length(Initial_position), ...
    Lower_bound, Upper_bound, options);

        
    function stop = myoutput(optimValues,state)
        stop = false;

    end
    
    function J = ESPM_fit_geom1(x) 
        x_opt_raw(1)=x(1);
        x_opt_raw(2)=x(2);
        
        x_opt_raw(3)=x(3);
        x_opt_raw(4)=x(4);
        
        x_opt_raw(9)=x(5);
        x_opt_raw(10)=x(6);
        
        x_opt_raw(11)=x(7);
        x_opt_raw(12)=x(8);
        
        x_opt_raw(15)=x(9);
        x_opt_raw(16)=x(10);
        
        x_opt_raw(17)=x(11);
        x_opt_raw(18)=x(12);
        
        x_opt_raw(20)=x(13);
        x_opt_raw(22)=x(14);
        
        x_opt_raw(23)=x(15);
        x_opt_raw(28)=x(16);
        x_opt_raw(25:27)=x(17:19);
        warning off
        try
      %%  Simulation
        [V_cell,soc_bulk_p,soc_bulk_n]=fun_espm(x_opt_raw,I_data,t_data,SOC_IC,profile_flag,int_method,model_type,battery_type);

           if sum(isreal(V_cell)) == 0 
               J = 10^3;
           else
               J=rms((V_cell-refV)./refV)+rms(refsocp-soc_bulk_p)+rms(refsocn-soc_bulk_n);
           end

        catch
            J = 10^3;
        end
            
    end
end




function [V_cell,soc_bulk_p,soc_bulk_n,param] = fun_espm(x_opt,I_data,t_data,SOC_IC,profile_flag,int_method,model_type,battery_type)






T_amb=23;



%% =====Run SPM or ESPM

[V_cell, R_l, T_core, T_surf,soc_bulk_n, soc_bulk_p, cs_n, cs_p,...
          V_oc,param,ocp_p,ocp_n,eta_p,eta_n,ce,eta_ele] = ESPM_main(x_opt,t_data,I_data,SOC_IC,T_amb,profile_flag,int_method,model_type,battery_type);

end



      
