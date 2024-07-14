% -------------------------------------------------------------------------
% Function for Up charge
% -------------------------------------------------------------------------
function y = U_p_cc_NMC_Pansonic(x)
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