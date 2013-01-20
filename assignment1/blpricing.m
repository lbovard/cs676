function op_price=blpricing(K,S0,r,tfinal,sigma,dt,option,gamma)
%%Binomial Lattice model for European options
%K ~ strike price
%S0 ~ initial price
%r ~ interest rate
%tfinal ~ final time (years)
%sigma ~ volatility
%dt ~ timestep size
%option ~ 0 for call, 1 for put
%gamma ~ power exponent

%size of lattice, rounded to nearest integer
N=round(tfinal/dt);
%Cox parameters
u=exp(sigma*sqrt(dt));
d=1/u;
%modified probability
qstar=(exp(r*dt)-d)/(u-d);
%precompute expectation values
pup=exp(-r*dt)*qstar;
pdown=exp(-r*dt)*(1-qstar);
%allocate V
V=zeros(N+1,1);
V(:)=S0*exp((2*(0:1:N)-N)*sigma*sqrt(dt));

if option==0
    %payoff for call
    V(:)=max(V(:)-ones(length(V(:)),1)*K,0).^gamma';
    %disp('Computing call option')
else
    %payoff for put
    V(:)=max(-V(:)+ones(length(V(:)),1)*K,0).^gamma';
    %disp('Computing put option')
end

%there is some optimisation in this loop since a lot of entries are 0 
%so calculations could be reduced. I played around with this and got huge
%performance increases but it was, for some bizarre reason, ill-condition
%wrt the initial interest rate. e.g. I changed r from 0.02 to 0.03 and the
%super optimised code went from ~ 8s to ~20s!
for n=[N:-1:1]
		V(1:n)=pup*V(2:(n+1))+pdown*V(1:n);
end
op_price=V(1);
end
