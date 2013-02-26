clear all;
close all;
K=100;
%call option
callpayoff=@(S) max(S-K,0);
putpayoff=@(S) max(K-S,0);
N=10;
[V0,L]=binomialDelta(100,0.01,0.2,1,N,callpayoff);

%alphat=binomialDelta(L(4).alpha,L(4).S,97);