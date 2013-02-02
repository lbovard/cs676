function Vout=mcbarrier(S0,Su,K,T,dt,r,sigma,M)
% Computes Monte-Carlo simulation for up-and-out barrier call options
% S0 ~ initial price; Su ~ upper barrier; K ~ strike price; T ~ final time
% dt ~ time-step; r ~ interest rate; sigma ~ volatility
% M ~ number of simulations

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
    S(1)=S0;
    %generate exact solution to risk neutral model
    F=exp((r-0.5*sigma^2)*dt+sigma*sqrt(dt)*randn(N+1,1));
    
    %generate stock price
    for i=2:N+1
        S(i)=S(i-1)*F(i);
    end
    %faster to check if barrier was hit outside loop due to vectorisation
    if(S(S>Su))
        %barrier was hit
        payoff(k)=0;
    else
        %barrier wasn't hit
        payoff(k)=max(S(end)-K,0);
    end
end 
Vout=exp(-r*T)*sum(payoff)/M;
end