%%Question 7 analysis
clear all;close all;
V=dlmread('q7_data.dat');
Vold=dlmread('q6b_data.dat');
%from part 6a
V_exact=0.401784554490646;
Vold=abs(Vold-V_exact)/V_exact;
%simulation data
dt=[5,1,0.1]/250;
M=[1000,2000,5000,10000,20000,50000,100000,500000,1000000];
%M=[1000,2000,5000,10000,20000,50000,100000];
figure
semilogx(M,V(1,:),'^g-')
hold on
semilogx(M,V(2,:),'*-')
semilogx(M,V(3,:),'ok-')
semilogx(M,V_exact*ones(length(M),1),'-r')
legend(num2str(dt(1)),num2str(dt(2)),num2str(dt(3)))
axis([0 M(end) 0.35 0.5])
xlabel('M')
ylabel('V_{0}')
title('Monte Carlo Simulation of Up-and-Out Option')

V=abs(V-V_exact)/V_exact;
figure
loglog(1./dt,Vold(:,end),'*-')
hold on
loglog(1./dt,V(:,end),'o-')
loglog(1./dt,0.1*dt,'-r')
loglog(1./dt,sqrt(dt),'--')
xlabel('N')
ylabel('Absolute Error')
title('Error in Monte Carlo Simulation of Up-and-Out Option')
legend('Old MC','New MC','Ref. Slope 1','Ref. Slope 1/2','location','SouthWest')

figure
%time data from simulations
time_old=[2.189552,6.236008,51.276373];
time_new=[5.446103,18.841457,147.148798];
loglog(1./dt,time_old,'*-')
hold on
loglog(1./dt,time_new,'o-')
xlabel('N')
ylabel('CPU Time')
title('Compute time of Monte Carlo Simulation of Up-and-Out Option')
legend('Old MC','New MC','location','NorthWest')