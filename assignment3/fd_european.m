function V=fd_european(sigma,r,T,K,S,ts_type,option_type,num_steps)
% Solve the Black-Scholes equation via finite difference
% sigma ~ volatility, r ~ interest rate, T ~final time, K ~ strike 
% S ~ asset grid
% ts_type ~ time stepping time (0=implicit,1=Crank-Nicolson,2=Rannacher)
% option_type ~ 0 call, 1 put
% num_steps ~ number of time steps 

%initial data
dt=T/num_steps;
m=length(S);

%set-up option price
if(option_type==0)
    V0=max(S-K,0)';
else
    V0=max(K-S,0)';
end

%% Set-up Matrix and Precondition
% define starting and final elements of grid, excluding boundary points
a=2;
b=m-1;
%set-up co-efficients. Note for simple BS we only need alpha_for,beta_for, alpha_cen,beta_cen
alpha_for=sigma^2*(S(a:b).^2)./((S(a:b)-S(a-1:b-1)).*(S(a+1:b+1)-S(a-1:b-1)));
beta_for=sigma^2*(S(a:b).^2)./((S(a+1:b+1)-S(a:b)).*(S(a+1:b+1)-S(a-1:b-1)))+r*S(a:b)./(S(a+1:b+1)-S(a:b));
alpha_back=sigma^2*(S(a:b).^2)./((S(a:b)-S(a-1:b-1)).*(S(a+1:b+1)-S(a-1:b-1))) - r*S(a:b)./(S(a:b)-S(a-1:b-1));
beta_back=sigma^2*(S(a:b).^2)./((S(a+1:b+1)-S(a:b)).*(S(a+1:b+1)-S(a-1:b-1)));
alpha_cen=sigma^2*(S(a:b).^2)./((S(a:b)-S(a-1:b-1)).*(S(a+1:b+1)-S(a-1:b-1))) - r*S(a:b)./(S(a+1:b+1)-S(a-1:b-1));
beta_cen=sigma^2*(S(a:b).^2)./((S(a+1:b+1)-S(a:b)).*(S(a+1:b+1)-S(a-1:b-1))) + r*S(a:b)./(S(a+1:b+1)-S(a-1:b-1));

%%precondition matrix elements to ensure alpha,beta >=0
f=alpha_for >=0 & beta_for >=0;
alpha_up(f)=alpha_for(f);
beta_up(f)=beta_for(f);
alpha_up(~f)=alpha_back(~f);
beta_up(~f)=beta_back(~f);
%get central co-efficients >=0
f= alpha_cen >=0 & beta_cen >=0 ;
alpha(f)=alpha_cen(f);
beta(f)=beta_cen(f);
%otherwise use forward co-efficients
alpha(~f)=alpha_up(~f);
beta(~f)=beta_up(~f);

%% TIME STEPPING 

if(ts_type==0)    
    
    %%IMPLICIT TIME STEPPING
    % now set-up M for MV^{n+1} = V^{n}
    M=sparse([2:m-1,1,m],[1:m-2,1,m],[-dt*alpha,0,0])+sparse([2:m-1,1,m],[3:m,1,m],[-dt*beta,0,0]) ... 
        + sparse(1:m,1:m,[1+r*dt,dt*alpha+dt*beta+r*dt+1,1]);

    %if put option value is 0 at boundary but that means det =0 
    % so must remove otherwise things blow up
    V=V0;
    if(option_type==1) 
        M=M(1:end-1,1:end-1);
        V=V0(1:end-1);
    end
    %do LU factorisation
    [L,U]=lu(M);
    for i=1:num_steps
        y=L\V;
        V=U\y;
    end

    %put in the other boundary, which is 0
    if(option_type==1)
        V(m)=0;
    end

elseif(ts_type==1)
    %%CRANK-NICOLSON TIME STEPPING
    %MUST SOLVE (1+M)V^{n+1}=(1-M)V^{n}

    %set up 1+M
    M=sparse([2:m-1,1,m],[1:m-2,1,m],[-dt*alpha/2,0,0])+sparse([2:m-1,1,m],[3:m,1,m],[-dt*beta/2,0,0]) ... 
        + sparse(1:m,1:m,[r*dt/2,dt*alpha/2+dt*beta/2+r*dt/2,1]);
    M=sparse(1:m,1:m,1)+M;

    %if put option must remove last entry to avoid errors
    if(option_type==1) 
        M=M(1:end-1,1:end-1);
    end
    %factorise 1+M 
    [L,U]=lu(M);

    %now define 1-M
    M=sparse([2:m-1,1,m],[1:m-2,1,m],[-dt*alpha/2,0,0])+sparse([2:m-1,1,m],[3:m,1,m],[-dt*beta/2,0,0]) ... 
        + sparse(1:m,1:m,[r*dt/2,dt*alpha/2+dt*beta/2+r*dt/2,1]);
    M=sparse(1:m,1:m,1)-M;

    V=V0(1:end-1);
    %again fix issue with last entry for put options 
    if(option_type==1) 
        M=M(1:end-1,1:end-1);
        V=V0(1:end-1);
    end
    for i=1:num_steps
        y=L\(M*V);
        V=U\y;
    end
    if(option_type==1)
        V(m)=0;
    end

else
    %%RANNACHER TIME STEPPING 


    %DO IMPLICIT TIME STEPS 
    %now set-up M for MV^{n+1} = V^{n}
    M=sparse([2:m-1,1,m],[1:m-2,1,m],[-dt*alpha,0,0])+sparse([2:m-1,1,m],[3:m,1,m],[-dt*beta,0,0]) ... 
        + sparse(1:m,1:m,[1+r*dt,dt*alpha+dt*beta+r*dt+1,1]);

    %if put option value is 0 at boundary but that means det =0 
    % so must remove otherwise things blow up
    V=V0;
    if(option_type==1) 
        M=M(1:end-1,1:end-1);
        V=V0(1:end-1);
    end
    %do LU factorisation
    [L,U]=lu(M);

    %do implicit time stepping 
    for i=1:2
        y=L\V;
        V=U\y;
    end
    %put in the other boundary, which is 0
    if(option_type==1)
        V(m)=0;
    end

    %now do CN 
    M=sparse([2:m-1,1,m],[1:m-2,1,m],[-dt*alpha/2,0,0])+sparse([2:m-1,1,m],[3:m,1,m],[-dt*beta/2,0,0]) ... 
        + sparse(1:m,1:m,[r*dt/2,dt*alpha/2+dt*beta/2+r*dt/2,1]);
    M=sparse(1:m,1:m,1)+M;

    %if put option must remove last entry to avoid errors
    if(option_type==1) 
        M=M(1:end-1,1:end-1);
    end
    %factorise 1+M 
    [L,U]=lu(M);

    %now define 1-M
    M=sparse([2:m-1,1,m],[1:m-2,1,m],[-dt*alpha/2,0,0])+sparse([2:m-1,1,m],[3:m,1,m],[-dt*beta/2,0,0]) ... 
        + sparse(1:m,1:m,[r*dt/2,dt*alpha/2+dt*beta/2+r*dt/2,1]);
    M=sparse(1:m,1:m,1)-M;

    %again fix issue with last entry for put options 
    if(option_type==1) 
        M=M(1:end-1,1:end-1);
        V=V(1:end-1);
    end

    %time-step CN
    for i=1:num_steps-2
        y=L\(M*V);
        V=U\y;
    end
    %add back in the endpoint for put options 
    if(option_type==1)
        V(m)=0;
    end

end


end