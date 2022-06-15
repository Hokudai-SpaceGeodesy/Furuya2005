function [u1,u2,u3]=qstherm(x1,x2,x3,xi1,xi2,xi3,t,Q)
% QSTHERM: compute quasi-static thermoelastic response due to an
% instantaneous heat source, Q.
%  >> [u1,u2,u3]=qstherm(x1,x2,x3,xi1,xi2,xi3,t,Q);
% Times steps "t" is given like, 
%  >> t=logspace(log10(10000),log10(365.25*24*3600*100),150);
% Initial heat source Q is assigned as "dT×dV"(K m^3); it is
% normalized by \rho * c. x3-axis is taken as positive downward.
% copyright (c) Masato Furuya, 2003
R=((x1-xi1)^2+(x2-xi2)^2+(x3-xi3)^2)^0.5;
R_=((x1-xi1)^2+(x2-xi2)^2+(-x3-xi3)^2)^0.5;
% material properties 
k=10^-5; % diffusivity !! Key parameter! 10^-6では遅過ぎる?
a=2*10^-5; % linear thermal expansivity
nu=0.25;    % Poisson ratio
m=a*(1+nu)/4/pi/(1-nu);
% multi-dimenstional "t" is allowed.
theta=4*k*t;
% commonly used function
f=erf(R./sqrt(theta))-2*R.*exp(-R^2./theta)./sqrt(pi*theta);
f_=erf(R_./sqrt(theta))-2*R_.*exp(-R_^2./theta)./sqrt(pi*theta);
% Infinity solutions
u1=((x1-xi1)/R^3)*f*m;
u2=((x2-xi2)/R^3)*f*m;
u3=((x3-xi3)/R^3)*f*m;
%
u1_=((x1-xi1)/R_^3)*f_*m;
u2_=((x2-xi2)/R_^3)*f_*m;
u3_=((-x3-xi3)/R_^3)*f_*m;
%
g_=2*R_*(1.+2*R_^2./3./theta).*exp(-R_^2./theta)./sqrt(pi*theta)- ...
   erf(R_./sqrt(theta));
e13_=(3*(x1-xi1)*(-x3-xi3)/R_^5)*g_*m;
e23_=(3*(x2-xi2)*(-x3-xi3)/R_^5)*g_*m;
e33_=(1/R_^3).*(1-3*(-x3-xi3)^2/R_^2).*f_*m+...
     m*4*(-x3-xi3)^2.*exp(-R_^2./theta)./R_^2./sqrt(pi)/theta.^1.5;
% Finally,...
u1=Q*(u1+(3-4*nu)*u1_-2*x3*e13_);
u2=Q*(u2+(3-4*nu)*u2_-2*x3*e23_);
u3=Q*(u3+(3-4*nu)*u3_+2*x3*e33_);