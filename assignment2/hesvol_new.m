function init = hesvol_new(K,S0,T,r,v_init,v_bar,eta,rho,lambda,M)
num_steps=100;
dt=T/100;
%cholesky factorisation of correlation matrix
Q=[1,rho;rho,1]; L=chol(Q); 

xold=log(S0)*ones(1,M);
vold=v_init*ones(1,M);
x=zeros(1,M);
v=zeros(1,M);

%number of Monte Carlo simulations
for i=2:num_steps+1
	%generate correlated random variables
	phi=L'*randn(2,M);
	%Heston stocastic volatility model 
    x=xold+r*dt-abs(vold/2)*dt+sqrt(abs(vold)*dt).*phi(1,:);
    v=(sqrt(abs(vold))+eta/2*sqrt(dt)*phi(2,:)).^2-lambda*(abs(vold)-v_bar)*dt-eta^2/4*dt;
    xold=x;
    vold=v;
		%x(1,i)=x(1,i-1)+r*dt-abs(x(2,i-1))/2*dt+sqrt(abs(x(2,i-1))*dt)*phi(1,i-1);
		%x(2,i)=(sqrt(abs(x(2,i-1)))+eta/2*sqrt(dt)*phi(2,i-1))^2-lambda*(abs(x(2,i-1))-v_bar)*dt-eta^2/4*dt;
end
Payoff=sum(max(exp(x)-K,0));
init=exp(-r*T)*Payoff/M;

end
