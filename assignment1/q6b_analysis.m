%%Question 6b analysis
clear all;close all;
V=dlmread('q6b_data.dat');
%from part 6a
V_exact=0.401784554490646;
%simulation data
dt=[5,1,0.1]/250;
M=[1000,2000,5000,10000,20000,50000,100000,200000,500000,1000000];
semilogx(M,V(1,:),'^g-')
hold on
semilogx(M,V(2,:),'*-')
semilogx(M,V(3,:),'ok-')
semilogx(M,V_exact*ones(length(M),1),'-r')
legend(num2str(dt(1)),num2str(dt(2)),num2str(dt(3)))
axis([0 M(end) 0.35 0.7])
xlabel('M')
ylabel('V_{0}')
title('Monte Carlo Simulation of Up-and-Out Option')