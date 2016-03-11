function PtW


r=log(2)
m=log(2)
phi=.05
beta_max=10
d=.1

Nmin=r*m/(r-d)\phi\beta_max;
Kmin=beta_max*phi/m*Nmin^2;
K=linspace(Kmin,100*Kmin,50);

N=sqrt(m/beta_max/phi*K);
V=(r-d)/phi-r/phi^2*m/beta_max./N;
DD_GrowthRate = r*(1-N./K);
Burst = beta_max*(1-N./K);
DeathRateByPhage=phi*V;

%figure(1)
%plot(N,Burst)
%title('PtW - Revolutionary')
%xlabel('Cells per ml')
%ylabel('Predicted burst size')
%grid on

%figure(2)
%plot(K,N)
%xlabel('Carrying Capacity')
%ylabel('Microbes')
%title('N vs K')
%grid on

%figure(3)
%plot(K,V)
%xlabel('Carrying Capacity')
%ylabel('VIRUSES')
%title('V vs K')
%grid on

%figure(4)
%plot(N./K,N)
%ylabel('N')
%xlabel('N/K')
%title('play')
%grid on

figure(5)
plot(N, DD_GrowthRate)
title('growth rate x microbes')
ylabel('DD_GrowthRate')
xlabel('Microbes')
grid on

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
%plot(V.*N, V./N)
%title('ER vs VMR')
%xlabel('ER')
%ylabel('VMR')
%grid on

%figure(9)
%plot(N, DeathRateByPhage)
%title('')
%xlabel('N')
%ylabel('DeathRateByPhage')
%grid on