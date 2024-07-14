% -------------------------------------------------------------------------
% Function for Un charge
% -------------------------------------------------------------------------
function y = U_n_cc_NMC_Pansonic(x)
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