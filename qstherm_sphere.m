function [u1,u2,u3]=qstherm_sphere(x1,x2,x3,xi1,xi2,xi3,t,Q,r)
% QSTHERM_SPHERE: compute quasi-static thermoelastic response due to an
% instantaneous heat source, Q in a spherical volume with radius r.
%  >> [u1,u2,u3]=qstherm_sphere(x1,x2,x3,xi1,xi2,xi3,t,Q,r);
% Times steps "t" is given like, 
%  >> t=logspace(log10(10000),log10(365.25*24*3600*100),150);
% Initial heat source Q is assigned as "dT x dV"(K m^3); it is
% normalized by \rho * c. x3-axis is taken as positive downward.<-- ??
% copyright (c) Masato Furuya, 2004-
R=((x1-xi1)^2+(x2-xi2)^2+(x3-xi3)^2)^0.5;
R_=((x1-xi1)^2+(x2-xi2)^2+(-x3-xi3)^2)^0.5;
% material properties 
k=10^-5; % diffusivity !! Key parameter! "10^-6" is too slow?
a=2*10^-5; % linear thermal expansivity
nu=0.25;    % Poisson ratio
m=a*(1+nu)/(1-nu);
% multi-dimenstional "t" is allowed.
theta=4*k*t;
% commonly used function 
sigma=3*Q/4/pi/r^3;
%
X=(R+r)./sqrt(theta);Y=(R-r)./sqrt(theta);
Phi_R=Q*m/4/pi/R^2+sigma*m*theta.^1.5/2/R^2.*(...
    ((X.^2/2-R*X./sqrt(theta)+R^2./theta).*ierfc(X)+(X.^3/6).*erfc(X)...
     -exp(-X.^2).*(X.^2+1)/6/sqrt(pi))-...
    ((Y.^2/2-R*Y./sqrt(theta)+R^2./theta).*ierfc(Y)+(Y.^3/6).*erfc(Y)...
     -exp(-Y.^2).*(Y.^2+1)/6/sqrt(pi)) ); 
Phi_RR=-2*sigma*m*r^3/3/R^3-sigma*m*theta.^1.5/R^3.*(...
    (0.5.*(X-R./sqrt(theta)).^2.*ierfc(X)+(X.^3/6-R^2.*X./2./theta+...
      R^3./2./theta.^1.5).*erfc(X)-((X.^2+1)./6).*exp(-X.^2)./sqrt(pi))-...
    (0.5.*(Y-R./sqrt(theta)).^2.*ierfc(Y)+(Y.^3/6-R^2.*Y./2./theta+...
      R^3./2./theta.^1.5).*erfc(Y)-((Y.^2+1)./6).*exp(-Y.^2)./sqrt(pi)));
% Infinity solutions 
u1=Phi_R*(x1-xi1)/R;
u2=Phi_R*(x2-xi2)/R;
u3=Phi_R*(x3-xi3)/R;
% conjugate derivatives
X=(R_+r)./sqrt(theta);Y=(R_-r)./sqrt(theta);
Phi_R_=Q*m/4/pi/R_^2+sigma*m*theta.^1.5/2/R_^2.*(...
    ((X.^2/2-R_*X./sqrt(theta)+R_^2./theta).*ierfc(X)+(X.^3/6).*erfc(X)...
     -exp(-X.^2).*(X.^2+1)/6/sqrt(pi))-...
    ((Y.^2/2-R_*Y./sqrt(theta)+R_^2./theta).*ierfc(Y)+(Y.^3/6).*erfc(Y)...
     -exp(-Y.^2).*(Y.^2+1)/6/sqrt(pi)) ); 
Phi_RR_=-2*sigma*m*r^3/3/R_^3-sigma*m*theta.^1.5/R_^3.*(...
    (0.5*(X-R_./sqrt(theta)).^2.*ierfc(X)+(X.^3/6-R_^2.*X./2./theta+...
      R_^3./2./theta.^1.5).*erfc(X)-((X.^2+1)./6).*exp(-X.^2)./sqrt(pi))-...
    (0.5*(Y-R_./sqrt(theta)).^2.*ierfc(Y)+(Y.^3/6-R_^2.*Y./2./theta+...
      R_^3./2./theta.^1.5).*erfc(Y)-((Y.^2+1)./6).*exp(-Y.^2)./sqrt(pi)));
u1_=Phi_R_*(x1-xi1)/R_;
u2_=Phi_R_*(x2-xi2)/R_;
u3_=Phi_R_*(-x3-xi3)/R_;
%
e13_=(Phi_RR_+Phi_R_/R_)*(x1-xi1)*(-x3-xi3)/R_^2;
e23_=(Phi_RR_+Phi_R_/R_)*(x2-xi2)*(-x3-xi3)/R_^2;
e33_=Phi_R_*(1-(-x3-xi3)^2/R_^2)/R_+Phi_RR_*(-x3-xi3)^2/R_^2;
% Finally,...
% In case of spherical volume, Q should not be multiplied; it was
% taken into account through sigma or Q.
if (R > r);
  u1=(u1+(3-4*nu)*u1_-2*x3*e13_);
  u2=(u2+(3-4*nu)*u2_-2*x3*e23_);
  u3=(u3+(3-4*nu)*u3_+2*x3*e33_);
else;
  u1=NaN;
  u2=NaN;
  u3=NaN;
end