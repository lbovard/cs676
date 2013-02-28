function [op_price,L]=binomialDelta(S0,r,sigma,T,N,fpayoff)
%%binomial tree hedging
dt=T/N;
%binomial model parameters
u=exp(sigma*sqrt(dt));d=1/u;qstar=(exp(r*dt)-d)/(u-d);
%precompute expectation values
pup=exp(-r*dt)*qstar;pdown=exp(-r*dt)*(1-qstar);
%allocate arrays
V=zeros(N+1,1); S=cell(1,N+1);alpha=S;
%initialise binomial stock prices
for n=1:N+1
    S{n}=S0*exp((2*((1:n)-1)-(n-1))*sigma*sqrt(dt));
end
%initialise option price
V(:)=S{N+1};
V=fpayoff(V');
%do binomial tree pricing and hedging
for n=N:-1:1
        alpha{n}=(V(2:(n+1))-V(1:n))./((u-d)*S{n});
		V(1:n)=pup*V(2:(n+1))+pdown*V(1:n);
end
op_price=V(1);
L=struct('S',S,'alpha',alpha);
end
