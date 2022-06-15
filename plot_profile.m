function plot_profile(func,t_in,unit,scale)

if ~exist("func", "var") func="sphere"; end
if ~exist("t_in", "var") t_in=1; end
if ~exist("unit", "var") unit="years"; end
if ~exist("scale", "var") scale=0; end

if (strcmp(unit,"months"));
  t=t_in*3600*24*12;
elseif (strcmp(unit,"years"));
  t=t_in*3600*24*365.25;
else (strcmp(unit,"seconds"));
  t=t_in;
end

x=[-99:2:100]; y=[-99:2:100]; z=[-198:2:0];
U1.d=zeros(10000,1);U2.d=zeros(10000,1);U3.d=zeros(10000,1);
depth=-50; q=5.*10^9; r=20;
p=9*10^10;

for i=1:100
  for j=1:100
    if (strcmp(func,"point"));
      [U1.d(j+100*(i-1),:),U2.d(j+100*(i-1),:),U3.d(j+100*(i-1),:)]=qstherm(0,y(i),z(j),0,0,depth,t,q);
    elseif (strcmp(func,"shell"));
      [U1.d(j+100*(i-1),:),U2.d(j+100*(i-1),:),U3.d(j+100*(i-1),:)]=qstherm_shell(0,y(i),z(j),0,0,depth,t,q,r);
    elseif (strcmp(func,"sphere"));
      [U1.d(j+100*(i-1),:),U2.d(j+100*(i-1),:),U3.d(j+100*(i-1),:)]=qstherm_sphere(0,y(i),z(j),0,0,depth,t,q,r);
    elseif (strcmp(func,"mogi"));
      [U1.d(j+100*(i-1),:),U2.d(j+100*(i-1),:),U3.d(j+100*(i-1),:)]=mogi(0,y(i),z(j),0,0,depth,t,p,r);
    end
  end
end

[U1a,U2a,U3a]=qs_snapshot(y,z,U1,U2,U3,1);
quiver(y,z,U2a,U3a,scale);
if (~strcmp(func,"mogi"));
  title([num2str(t_in,"%.1f") " " unit]);
else
  title("Mogi source");
end
xlim([-100 100])
ylim([-200 50])
axis equal;