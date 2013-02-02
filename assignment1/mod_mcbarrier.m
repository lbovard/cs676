function Vout=mod_mcbarrier(S0,Su,K,T,dt,r,sigma,M)
% Computes Monte-Carlo simulation for up-and-out barrier call options
% S0 ~ initial price; Su ~ upper barrier; K ~ strike price; T ~ final time
% dt ~ time-step; r ~ interest rate; sigma ~ volatility
% M ~ number of simulations
%testing data
%S0=100;K=95;Su=110;r=0.01;sigma=0.2;T=1;
%dt=0.1/250;
%M=5000;
%initial data
N=round(T/dt);
%number of simulations

%set-up stuff
payoff=zeros(1,M);

%Monte Carlo simulation
%let's parallelise this because MATLAB makes it so easy!
%plus it makes code zoom zoom
parfor k=1:M
    S=zeros(N+1,1);
    P=S;
    S(1)=S0;
    %generate exact solution to risk neutral model
    F=exp((r-0.5*sigma^2)*dt+sigma*sqrt(dt)*randn(N+1,1));
    %generate stock price
    for i=2:N+1
        S(i)=S(i-1)*F(i);
        %brownian bridge probability
        P(i)=exp(-2*(Su-S(i-1))*(Su-S(i))/(sigma^2*S(i-1)^2*dt));
    end
    %faster to check if barrier was hit outside loop due to vectorisation
    %define array to check if Brownian bridge probablity occurs
    u=rand(N+1,1);
    %check to see if the barrier price is reached
    A=S(S>Su);
    %check to see if the bridge was crossed
    B=u(P>u);
    if(~isempty(A) || ~isempty(B))
        %barrier was hit
        payoff(k)=0;
    else
        %barrier wasn't hit
        payoff(k)=max(S(end)-K,0);
    end
end 
Vout=exp(-r*T)*sum(payoff)/M;
end