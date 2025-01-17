%%Question 6a 
% Evaluating barrier options using Monte Carlo simulation
% with exact solution from Black-Scholes equation for 
% up-and-out European options
clear all;close all;
%set some variables
T=1;t=0;sigma=0.2;r=0.01;K=95;Su=110;
S=80:2:110;
%norm-cdf parameters
d1=(log(S/K)+(r+0.5*sigma^2)*(T-t))/(sigma*sqrt(T-t));
d2=d1-sigma*sqrt(T-t);
d3=(log(S/Su)+(r+0.5*sigma^2)*(T-t))/(sigma*sqrt(T-t));
d4=d3-sigma*sqrt(T-t);
d5=(log(S/Su)-(r-0.5*sigma^2)*(T-t))/(sigma*sqrt(T-t));
d6=d5-sigma*sqrt(T-t);
d7=(log(S*K/Su^2)-(r-0.5*sigma^2)*(T-t))/(sigma*sqrt(T-t));
d8=d7-sigma*sqrt(T-t);

%compute exact initial price
Vout=S.*(normcdf(d1)-normcdf(d3)-(Su./S).^(1+2*r/sigma^2).*(normcdf(d6)-normcdf(d8))) - ...
    K*exp(-r*(T-t))*(normcdf(d2)-normcdf(d4)-(Su./S).^(-1+2*r/sigma^2).*(normcdf(d5)-normcdf(d7)));

plot(S,Vout,'*-')
title('European up-and-out barrier call option')
xlabel('Initial Price')
ylabel('V_{out}(S,0)')
axis([80 110 0 0.6])