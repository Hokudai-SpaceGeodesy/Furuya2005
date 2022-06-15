t=logspace(log10(10000),log10(365.25*24*3600*1000),150);
t_yr=t/365.25/24/3600;

x1=0; x2=0; x3=0;
xi1=0; xi2=0; xi3=-50;
Q=5.*10^9;
r=20;

[u1,u2,u3a_0]=qstherm(x1,x2,x3,xi1,xi2,xi3,t,Q);
[u1,u2,u3b_0]=qstherm_shell(x1,x2,x3,xi1,xi2,xi3,t,Q,r);
[u1,u2,u3c_0]=qstherm_sphere(x1,x2,x3,xi1,xi2,xi3,t,Q,r);

x3=-20;
[u1,u2,u3a_20]=qstherm(x1,x2,x3,xi1,xi2,xi3,t,Q);
[u1,u2,u3b_20]=qstherm_shell(x1,x2,x3,xi1,xi2,xi3,t,Q,r);
[u1,u2,u3c_20]=qstherm_sphere(x1,x2,x3,xi1,xi2,xi3,t,Q,r);

subplot(2,1,1)
hold on;
title("Vertical displacement at the surface");
semilogx(t_yr, u3a_0, '-', 'Color', [0,0,0]);
semilogx(t_yr, u3b_0, '-.', 'Color', [0,0,0]);
semilogx(t_yr, u3c_0, ':', 'Color', [0,0,0]);
xlim([10^-3 10^3])
ylim([0 20])
xlabel("Time (years)")
ylabel("Displacement (m)")
legend("Point","Shell","Volume")
hold off;

subplot(2,1,2)
hold on;
title("Vertical displacement at depth of 20 m");
semilogx(t_yr, u3a_20, '-', 'Color', [0,0,0]);
semilogx(t_yr, u3b_20, '-.', 'Color', [0,0,0]);
semilogx(t_yr, u3c_20, ':', 'Color', [0,0,0]);
xlim([10^-3 10^3])
ylim([0 30])
xlabel("Time (years)")
ylabel("Displacement (m)")
legend("Point","Shell","Volume")
hold off;