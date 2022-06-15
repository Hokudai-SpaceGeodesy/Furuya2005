function [u1,u2,u3]=qstherm_shell(x1,x2,x3,xi1,xi2,xi3,t,Q,r)
% QSTHERM_SHELL: compute quasi-static thermoelastic response due to an
% instantaneous heat source, Q in a spherical surface with radius r.
%  >> [u1,u2,u3]=qstherm_shell(x1,x2,x3,xi1,xi2,xi3,t,Q,r);
% Times steps "t" is given like, 
%  >> t=logspace(log10(10000),log10(365.25*24*3600*100),150);
% Initial heat source Q is assigned as "dT x dV"(K m^3); it is
% normalized by \rho * c. x3-axis is taken as positive downward.
% copyright (c) Masato Furuya, 2004-
R=((x1-xi1)^2+(x2-xi2)^2+(x3-xi3)^2)^0.5;
R_=((x1-xi1)^2+(x2-xi2)^2+(-x3-xi3)^2)^0.5;
% material properties 
k=10^-5; % diffusivity !! Key parameter! "10^-6" is too slow?
a=2*10^-5; % linear thermal expansivity
nu=0.25;    % Poisson ratio
m=a*(1+nu)/(1-nu); % <- correct definition for m in this case.
                   %    Do not confuse with that in point source.
% multi-dimenstional "t" is allowed.
theta=4*k*t;
% commonly used function 
sigma=Q/4/pi/r^2;
%
X=(R+r)./sqrt(theta);Y=(R-r)./sqrt(theta);
Phi_R=1+(1/2/r)*sqrt(theta).*(...
    ierfc(X)+R*erfc(X)./sqrt(theta)-...
   (ierfc(Y)+R*erfc(Y)./sqrt(theta)));
Phi_R=Phi_R*sigma*r^2*m/R^2;
Phi_RR=1+(1/2/r)*sqrt(theta).*(...
 (ierfc(X)+R*erfc(X)./sqrt(theta)+(R./sqrt(pi)./theta.^1.5).*exp(-X.^2))-...
 (ierfc(Y)+R*erfc(Y)./sqrt(theta)+(R./sqrt(pi)./theta.^1.5).*exp(-Y.^2)));
Phi_RR=-Phi_RR*2*sigma*m*r^2/R^3;
% Infinity solutions 
u1=Phi_R*(x1-xi1)/R;
u2=Phi_R*(x2-xi2)/R;
u3=Phi_R*(x3-xi3)/R;
% conjugate derivatives
X=(R_+r)./sqrt(theta);Y=(R_-r)./sqrt(theta);
Phi_R_=1+(1/2/r)*sqrt(theta).*(...
    ierfc(X)+R_*erfc(X)./sqrt(theta)-...
   (ierfc(Y)+R_*erfc(Y)./sqrt(theta)));
Phi_R_=Phi_R_*sigma*r^2*m/R_^2;
Phi_RR_=1+(1/2/r)*sqrt(theta).*(...
 (ierfc(X)+R*erfc(X)./sqrt(theta)+(R_./sqrt(pi)./theta.^1.5).*exp(-X.^2))-...
 (ierfc(Y)+R*erfc(Y)./sqrt(theta)+(R_./sqrt(pi)./theta.^1.5).*exp(-Y.^2)));
Phi_RR_=-Phi_RR_*2*sigma*m*r^2/R_^3;
u1_=Phi_R_*(x1-xi1)/R_;
u2_=Phi_R_*(x2-xi2)/R_;
u3_=Phi_R_*(-x3-xi3)/R_;
%
e13_=(Phi_RR_+Phi_R_/R)*(x1-xi1)*(-x3-xi3)/R_^2;
e23_=(Phi_RR_+Phi_R_/R)*(x2-xi2)*(-x3-xi3)/R_^2;
e33_=Phi_R_*(1-(-x3-xi3)^2/R_^2)/R_+Phi_RR_*(-x3-xi3)^2/R_^2;
% Finally,...
% In case of spherical volume, Q should not be multiplied; it is
% already taken into account through sigma or Q.
if (R > r);
  u1=(u1+(3-4*nu)*u1_-2*x3*e13_);
  u2=(u2+(3-4*nu)*u2_-2*x3*e23_);
  u3=(u3+(3-4*nu)*u3_+2*x3*e33_);
else;
  u1=NaN;
  u2=NaN;
  u3=NaN;
end