%%Q5 properties of Brownian Motion
clear all;close all;
N=[100,1000];
figure
for i=1:2
    h(i)=subplot(1,2,i);
    dt=1/N(i);
    x=0:dt:1;
    for j=1:10
        %generate random numbers
        z=sqrt(dt)*randn(N(i)+1,1);
        %do cumulative sum
        z=cumsum(z);
        plot(x,z,'k')
        hold on
    end
    axis square
    %axis([0 1 -3 3])
    xlabel('Time')
    ylabel('Z_{t}')
    title(strcat('Brownian Motion for N=',num2str(N(i))))
end
linkaxes([h(1),h(2)],'y')

figure
N=[100,1000];
for i=1:2
    h(i)=subplot(1,2,i);
    dt=1/N(i);
    x=0:dt:1;
    for j=1:10
        %do new modified Brownian motion
        z=zeros(N(i)+1,1);
        z(1)=sqrt(dt*randn^2);
        l=dt*randn(N(i)+1,1).^2;
        for k=2:(N(i)+1)
            z(k)=sqrt(z(k-1)^2+l(k));
        end
        plot(x,z,'k')
        hold on
    end
    axis square
    xlabel('Time')
    ylabel('Z_{t}')
    title(strcat('Modified Brownian Motion for N=',num2str(N(i))))
end
