%%Binomial Lattice model for pricing
close all;
clear all;
%time step
dt=1/3;
%final time (years)
tfinal=1;
%number of lattices
N=tfinal/dt;
%strike price
K=100;
%initial price
S0=100;
%volatility
sigma=0.2;
%interest
r=0.01;

%Cox parameters
u=exp(sigma*sqrt(dt));
d=1/u;

%modified probability
qstar=(exp(r*dt)-d)/(u-d);

%allocate V
V=zeros(N+1,N+1);
S=zeros(N+1,N+1);
%final value of asset
S(:,N+1)=S0*exp((2*(0:1:N)-N)*sigma*sqrt(dt));
%payoff 
V(:,N+1)=max(S(:,N+1)-ones(length(S(:,N+1)),1)*K,0)';

%do the looping backwards
for n=[N:-1:1]
	for jj=[0:1:n-1] 
		j=jj+1;
		S(j,n)=exp(-r*dt)*(qstar*S(j+1,n+1)+(1-qstar)*S(j,n+1));		
%		V(j,n)=exp(-r*dt)*(qstar*V(j+1,n+1)+(1-qstar)*V(j,n+1));		
	end
end
V(1,1)
