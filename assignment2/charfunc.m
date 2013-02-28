function y=charfunc(phi,lambda,vbar,rho,eta,r,tau,x,v,K,term)
    %phi is the dummy integration variable
    %lambda,vbar,rho,eta are parameters
    %term selections 1 or 2 characteristic function
    a=lambda*vbar;
    if(term==1)
        b=lambda-rho*eta;
        u=0.5;
    else
        b=lambda;
        u=-0.5;
    end
    %equation 17-18 of Heston 1993
    d=sqrt((rho*eta*phi*1i-b).^2-eta^2*(2*u*phi*1i-phi.^2));
    g=(b-rho*eta*phi*1i+d)./(b-rho*eta*phi*1i-d);
    C=r*phi*1i*tau+a/eta^2*(tau*(b-rho*eta*phi*1i-d)-2*log((exp(-d*tau)-g)./(1-g)));
    D=(b-rho*eta*phi*1i+d)/eta^2 .*((1-exp(d*tau))./(1-g.*exp(d*tau)));
    f=exp(C+D*v+1i*phi*x);
    y=real(exp(-1i*phi*log(K)).*f./(1i*phi));
end
