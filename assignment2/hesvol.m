function init = hesvol(K,S0,T,r,v_init,v_bar,eta,rho,lambda,M)
num_steps=100;
dt=T/100;
%cholesky factorisation of correlation matrix
Q=[1,rho;rho,1]; L=chol(Q); 

x=zeros(2,num_steps+1);
Payoff=0;
for m=1:M
	%generate correlated random variables
	phi=L'*randn(2,100);
	%initialise simulation
	x(1,1)=log(S0); x(2,1)=v_init;
	%Heston stocastic volatility model 
	for i=2:num_steps+1
		x(1,i)=x(1,i-1)+r*dt-abs(x(2,i-1))/2*dt+sqrt(abs(x(2,i-1))*dt)*phi(1,i-1);
		x(2,i)=(sqrt(abs(x(2,i-1)))+eta/2*sqrt(dt)*phi(2,i-1))^2-lambda*(abs(x(2,i-1))-v_bar)*dt-eta^2/4*dt;
	end
	Payoff=Payoff+max(exp(x(1,end))-K,0);
end
init=exp(-r*T)*Payoff/M;
%t=linspace(0,0.5,num_steps+1);
%plot(t,exp(x(1,:)))
%plot(t,x(2,:))

end
