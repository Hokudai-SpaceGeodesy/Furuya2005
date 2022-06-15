function [U1a,U2a,U3a]=qs_snapshot(x,y,U1,U2,U3,n)
% QS_SNAPSHOT: retrieve snapshot of temporal thermoelastic
% response, which has been computed like the following:
% >> t=logspace(log10(10000),log10(365.25*24*3600*100),150);
% >> x=[-99:2:100]; y=[-99:2:100]; z=[-198:2:0]
% >> U1.d=zeros(10000,150);U2.d=zeros(10000,150);U3.d=zeros(10000,150);
% >> for i=1:100
%       for j=1:100
%         [U1.d(j+100*(i-1),:),U2.d(j+100*(i-1),:),U3.d(j+100*(i-1),:)]=qstherm(x(i),y(j),0,0,0,-50,t,5*10^9);
%       end
%     end
% Usage:
% >> [U1a,U2a,U3a]=qs_snapshot(x,y,U1,U2,U3,n)
% >> quiver(x,y,U1a,U2a,0)
% >> axis equal
% where "n" stands for the index number of "t".
% copyright (c) Masato Furuya, 2004 - 
Nx=length(x);Ny=length(y);
U1a=zeros(Ny,Nx);U2a=zeros(Ny,Nx);U3a=zeros(Ny,Nx);
U=U1.d(:,n);V=U2.d(:,n);W=U3.d(:,n);
for i=1:Nx
  for j=1:Ny
    U1a(j,i)=U(j+Nx*(i-1),1);
    U2a(j,i)=V(j+Nx*(i-1),1);
    U3a(j,i)=W(j+Nx*(i-1),1);
  end
end