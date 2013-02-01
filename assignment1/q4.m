%%Run binomial lattice simulation for various dt
%for both call and put
clear all;
close all;
%initial parameters
K=100;S0=100;r=0.01;sigma=0.2;tfinal=1;
%number of times to sample 
N=6;
timesteps=0.01./2.^(0:(N-1));
%allocate arrays
call_data=zeros(4,length(timesteps));
call_data(1,:)=timesteps;
put_data=call_data;

%%Part A
%run simulation
i=1;
for dt=timesteps 
    call_data(2,i)=blpricing(K,S0,r,tfinal,sigma,dt,0,1);
    put_data(2,i)=blpricing(K,S0,r,tfinal,sigma,dt,1,1);
    i=i+1;
end

%compute change and ratio
for i=2:N
    call_data(3,i)=call_data(2,i)-call_data(2,i-1);
    put_data(3,i)=put_data(2,i)-put_data(2,i-1);
end
for i=3:N
    call_data(4,i)=(call_data(2,i-1)-call_data(2,i-2))/(call_data(2,i)-call_data(2,i-1));
    put_data(4,i)=(put_data(2,i-1)-put_data(2,i-2))/(put_data(2,i)-put_data(2,i-1));
end

%%Part B
strikes=60:10:100;
dt=0.0025;
gamma=2 ;
i=1;
for K=strikes
    power_put(1,i)=blpricing(K,S0,r,tfinal,sigma,dt,1,gamma);
    power_put(2,i)=blpricing(K,S0,r,tfinal,sigma,dt,1,1);
    i=i+1;
end
plot(strikes,power_put(1,:),'*-')
xlabel('Strike Price')
ylabel('Initial Put Price')
hold on
plot(strikes,power_put(2,:),'or-')
title('European Power Put')
legend('\gamma=2','\gamma=1','Location','NorthWest')