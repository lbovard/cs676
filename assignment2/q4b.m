close all;
clear all;
K=100;S0=100;T=0.5;r=0.02;vinit=0.0174;vbar=0.0354;eta=0.3877;rho=0.7165;lambda=1.3253;M=10000;

K=linspace(0.7*S0,1.2*S0,20);
init=hesvol(K,S0,T,r,vinit,vbar,eta,rho,lambda,M);

impv=blsimpv(S0,K,r,T,init);
figure
plot(K,impv,'*-')
xlabel('Strike K')
ylabel('\sigma implied')
grid on

K=100;
T=linspace(0.2,0.5,10);
init=zeros(1,length(T));impv=init;
for i=1:length(T)
    init(i)=hesvol(K,S0,T(i),r,vinit,vbar,eta,rho,lambda,M);
    impv(i)=blsimpv(S0,K,r,T(i),init(i));
end
figure
plot(T,impv,'*-')
xlabel('Strike K')
ylabel('\sigma implied')
grid on