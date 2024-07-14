% -------------------------------------------------------------------------
% Function for Up discharge
% -------------------------------------------------------------------------
function y = U_p_dis_NMC_Pansonic(x)
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