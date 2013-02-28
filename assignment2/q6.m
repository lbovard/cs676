clear all;
close all;tic
%compute the exact solution to the Heston stochastic volatility model
lambda=1.3253;vbar=0.0354;rho=0.7165;eta=0.3877;r=0.02;
tau=0.5;
K=100;
S0=100;
x=log(S0);
v=0.0174;

%define the characteristic functions
f1=@(xx)charfunc(xx,lambda,vbar,rho,eta,r,tau,x,v,K,1);
f2=@(xx)charfunc(xx,lambda,vbar,rho,eta,r,tau,x,v,K,2);

%integrate the characteristic functions from 0 to inf
%since f1,f2 formula have a 1/phi term matlab will divide by zero however
%the function is easily seen to be finite (via limits) so instead of doing
%some fancy tricks, just start from a number slightly above 0. The upper
%bound of inf can be handled by the formula above gives NaNs at some point
%so instead just set inf = 5000 which is pretty good because
%characteristic functions decay sufficiently quickly 
P1=quadgk(f1,1e-10,5000,'RelTol',1e-12,'AbsTol',1e-12)
P2=quadgk(f2,1e-10,5000,'RelTol',1e-12,'AbsTol',1e-12);

P1=0.5+P1/pi;
P2=0.5+P2/pi;

%exact formula
C=S0*P1-K*exp(-r*tau)*P2
toc
%do MC simulation
K=100;S0=100;T=0.5;r=0.02;vinit=0.0174;vbar=0.0354;eta=0.3877;rho=0.7165;lambda=1.3253;M=50000;
init=hesvol(K,S0,T,r,vinit,vbar,eta,rho,lambda,M)

%compare with heston numerical model
norm(C-init)

