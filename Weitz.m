function Weitz

a=.1
r=log(2)
m=log(2)
phi=.05
beta=10
d=.1


Kmin=1.8
% Kmin is chosen to insure V > 0. If you change the parameters,
% use fig 3 to adjust Kmin.
Kmax=10*Kmin

K=linspace(Kmin,Kmax,50);
N=(K/2/a).*(1-sqrt(1-4*m*a/phi/beta./K));
slope=-(a*d+(1-a)*r)*beta/a/m;
V=r/a/phi+slope*N;
DD_GrowthRate = r*(1-N./K);
DD_DeathRateByPhage=DD_GrowthRate-d;

figure(1)
plot(N,DD_GrowthRate)
title('Weitz and Dushoff')
xlabel('Cells per ml')
ylabel('Density-Dependent Growth Rate')
grid on

%figure(2)
%plot(K,N)
%xlabel('Carrying Capacity')
%ylabel('Microbes')
%title('N vs K')
%grid on

%figure(3)
%plot(N,DD_GrowthRate)
%xlabel('Microbes')
%ylabel('Density-dependent GrowthRate')
%title('')
%grid on

%figure(4)
%plot(N./K, DD_DeathRateByPhage)
%ylabel('death')
%xlabel('N/K')
%title('')
%grid on

%figure(5)
%plot(K,V./N)
%title('VMR vs K')
%xlabel('K')
%ylabel('VMR')
%grid on

figure(6)
plot(N,V)
title('V vs N')
xlabel('Microbes')
ylabel('Viruses')
grid on

%figure(7)
%plot(K,V.*N)
%title('ER vs K')
%xlabel('Carrying Capacity')
%ylabel('Encounter Rate (N*V)')
%grid on

%figure(8)
%plot(N, DD_GrowthRate)
%title('')
%ylabel('Density-dependent GrowthRate')
%xlabel('Microbes')
%grid on

%figure(9)
%plot(N, DD_DeathRateByPhage)
%title('')
%ylabel('Density-dependent DeathRate by Phage')
%xlabel('Microbes')
%grid on