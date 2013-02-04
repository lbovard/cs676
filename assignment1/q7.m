%%Question 7
%Monte Carlo simulation of up and out barrier option modified by Moon
close all;clear all;
%initial data
S0=100;K=95;Su=110;r=0.01;sigma=0.2;T=1;

%let's choose some dt, M
dt=[5,1,0.1]/250;
M=[1000,2000,5000,10000,20000,50000,100000,500000,1000000];

V=zeros(length(dt),length(M));
j=1;
for ts=dt
    i=1;
    disp(strcat('Computing timestep=',num2str(ts)))
    for mv=M
        disp(strcat('Computing simulation M=',num2str(mv)))
        tic
        V(j,i)=mod_mcbarrier(S0,Su,K,T,ts,r,sigma,mv);
        toc
        i=i+1;
    end
    j=j+1;
end

%save the data
dlmwrite('q7_data.dat',V,'precision',15);