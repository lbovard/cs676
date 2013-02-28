function [var,cvar]=dVaRCVaR(L,beta)
    L=sort(L);
    M=length(L);
    %since ib/M >= beta > (ib-1)/M  
    %implies beta*M <= ib && ib < ib*beta +1 
    %so ib = ceil(M*beta)
    ib=ceil(M*beta);
    var=L(ib);
    cvar=1/(1-beta)*((ib/M-beta)*L(ib)+1/M*sum(L(ib+1:M)));
end