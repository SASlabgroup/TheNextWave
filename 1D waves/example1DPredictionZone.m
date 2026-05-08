clear; close all; clc; 

% Normalize by wave length/period or something to give context to x,t


T = 17; % data collection duration
xi = 6; % buoy position
cg_min = 5; % minimum group velocity
cg_max = 10; % maximum group velocity

beta = cg_min/(cg_max - cg_min);
Ti = -beta*T;
Tf = T + beta*T;

t1 = linspace(Ti,0,100);
t2 = linspace(0,T,100);
t3 = linspace(T,Tf,100);

f1 = figure();
hold on
grid on
% for Ti < t < 0
plot((t1-T)*cg_min + xi,t1,'r')
plot(t1*cg_max + xi,t1,'b')

% for 0 < t < T
plot((t2-T)*cg_min + xi,t2,'r')
plot(t2*cg_min + xi,t2,'b')

% for T < t < Tf
plot((t3-T)*cg_max + xi,t3,'r')
plot(t3*cg_min + xi,t3,'b')

xline(xi,'k--')
xline(xi + beta*cg_max*T,'b--')
xline(xi - beta*cg_max*T,'r--')
yline(T,'k--')
yline(0,'k--')
yline(Ti,'b--')
yline(Tf,'r--')

% xlim([-1 1]*T*cg_max)
% ylim([Ti Tf])

xlabel('x [m]')
ylabel('t [s]')

title(strcat("T = ",num2str(T)," s, $\xi$ = ",num2str(xi)," m, $c_{g,min}$ = ",num2str(cg_min)," m/s, $c_{g,max}$ = ",num2str(cg_max)," m/s"))