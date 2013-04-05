clear all;
close all;

%initial data
sigma=0.38;r=0.025;T=0.5;K=100;


S=[0:0.1*K:0.4*K,0.45*K:0.05*K:0.8*K,0.82*K:0.02*K:0.9*K,0.91*K:0.01*K:1.1*K,1.12*K:0.02*K:1.2*K,1.25*K:.05*K:1.6*K,1.7*K:0.1*K:2*K,2.2*K,2.4*K,2.8*K,3.6*K,5*K,7.5*K,10*K];

time_step=2;
option_type=1;
num_step=25;
for k=1:9
    % exact solution
    [call,put]=blsprice(S,K,r,T,sigma);
    V=fd_european(sigma,r,T,K,S,time_step,option_type,num_step);
    L(k)=norm(put-V',inf);
    
    op=S;
    
    for i=1:length(S)-1
        P(2*i-1)=S(i);
        P(2*i)=(S(i)+S(i+1))/2;
    end
    S=P;
    num_step=num_step*2;
end
ratios=(L(1:(k-2))-L(2:(k-1)))./(L(2:(k-1))-L(3:k));

a=find(op==50);
b=find(op==150);

A=op(a:b);B=V(a:b);

figure
plot(A,B,'-')untitled.png
title('Put Option Value w/ Rannacher Timestepping')
ylabel('Value')
xlabel('Asset Price')
grid on

del=(B(2:end)-B(1:end-1))'./(A(2:end)-A(1:end-1));
figure
plot(A(1:end-1),del,'-')
title('Delta of Put Option w/ Rannacher Timestepping')
ylabel('Delta')
xlabel('Asset Price')
grid on
gamma=2*(del(2:end)-del(1:end-1))./(A(3:end)-A(2:end-1));

figure
plot(A(3:end),gamma,'-')
title('Gamma of Put Option w/ Rannacher Timestepping')
ylabel('Gamma')
xlabel('Asset Price')
grid on