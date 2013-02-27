clear all;
close all;
sigma=.3;mu=0.15;T=1.0;S0=100;r=0.02;K=90;
N=12;dt=T/N;
%number of MC simulations
M=10000;
%generate binomial model
payoff=@(S) max(S-K,0);
[V0,L]=binomialDelta(S0,r,sigma,T,N,payoff);
PL=zeros(1,M);
parfor k=1:M
    %initialise stock prices
    S=ones(1,N+1);S(1)=S0;
    S(2:end)=exp((mu-0.5*sigma^2)*dt+sigma*randn(1,N)*sqrt(dt));
    %compute recursive relationship
    S=cumprod(S);

    %initialise variables
    alphat=zeros(1,N);bank=alphat;A=exp(r*dt);
    %get first hedge position
    alphat(1)=L(1).alpha;

    %get hedging positions
    %couldn't vectorise this as interp1 complains about the arrays being
    %different sizes. Could be fixed with cleverness?
    %alphat contains hedging times at t_{0},t_{1},t_{2},...,t_{N-1}
    for i=2:N
        alphat(i)=interpDelta(L(i).alpha,L(i).S,S(i));
    end
    %initial bank
    bank(1)=V0-alphat(1)*S(1);
    %some pre-computing
    f=(alphat(2:N)-alphat(1:N-1)).*S(2:N);
    %rebalance bank
    %don't think this can be vectorised easily despite inhomogeneous recurrence
    %having an exact formula. But it's far faster just to compute the
    %recurrence 
    for i=2:N
        %bank(n)=A.^(n-1)*bank(1)-sum(A.^(n-1-(1:n-1)).*f(1:n-1));
       bank(i)=A*bank(i-1)-f(i-1);
    end
    %note that S and alphat/bank are different sizes
    %S(end) = stock price at t_{N}
    %alpha/bank(end) = hedging/bank at t_{N-1}
    port=-payoff(S(end))+alphat(end)*S(end)+bank(end)*A;

    PL(k)=exp(-r*T)*port/V0;
    
end