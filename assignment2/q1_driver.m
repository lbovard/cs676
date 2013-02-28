clear all;
close all;
sigma=0.3;r=0.02;T=1.0;K=90;S0=100;N=100;

callpayoff=@(S) max(S-K,0);
putpayoff=@(S) max(K-S,0);
[V0,L]=binomialDelta(S0,r,sigma,T,N,callpayoff);
St=linspace(80,140,100);
for i=1:length(St)
    alphat(i)=interpDelta(L(end-1).alpha,L(end-1).S,St(i));
end
alpha_exact=blsdelta(St,K,r,T/N,sigma);

plot(St,alpha_exact,'r*--','LineWidth',1,'MarkerSize',9)
hold on
plot(St,alphat,'^--','LineWidth',1,'MarkerSize',9)
xlabel('Asset price','FontSize',14)
ylabel('\alpha','FontSize',14)
title('Comparison of Exact and Binomial Tree Delta Hedging','FontSize',14)
axis([80 110 -0.05 1.05])
grid on
legend('Exact','Binomial','Location','NorthWest')
set(gca,'FontSize',14)
