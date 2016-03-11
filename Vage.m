function Vage

delta0=1;
D=.1;
mu0max=1;
m0=20;
beta0=1.5e-10;

sigma=.5;
rho=.9;
alpha=.2;
nu=0.9;

SRmax=400
numSRvalues=200;
% ******************************************
H0=(delta0+D)/((m0-1)*beta0);
Q=5e-8;
SRmin1=H0*Q
SRmin3=nu*mu0max*D/alpha/(nu*mu0max-D);
SRmin=max(SRmin1,SRmin3);
SRvalues=linspace(SRmin,SRmax,numSRvalues);
for k=1:length(SRvalues)
    SR=SRvalues(k);
    % ***** DILUTION LIMITATION *****
    factor=1-D/alpha/SR;
    n_d=floor(log(D/mu0max/factor)/log(nu));
    % *******************************
    % ***** CARRYING CAPACITY *****
    HT=SR/Q;
    CC(k)=HT;
    % ***** RESOURCE LIMITATION *****
    rho=rho+eps;
    aa=(1-rho)/(1-rho*sigma)*(HT/H0-1);
    n_r=floor(1-log(1+aa)/log(rho));
    n=n_r;
    if n>n_d
        n=n_d;
    end
    num_strains(k)=n;
    S=SR;
    mumax=mu0max*nu.^[0:n-1]';
    mu=mumax*alpha*S./(alpha*S+mumax);
    c=mu-D;
    % ADSORPTION MATRIX
    b=beta0*beta(rho,sigma,n);
    H=(delta0+D)/(m0-1)*((b')^(-1))*ones(n,1);
    
    
    V=(b^(-1))*c;
    H_total=sum(H);
    % SHANNON DIVERSITY
    p=H/H_total;
    ShannonDiversity(k)=-sum(p.*log(p));
    V_total=sum(V);
    VHR(k)=V_total/H_total;
    Hx(k)=H_total;
    Vx(k)=V_total;
    % AVERAGE INTRINSIC GROWTH RATE
    AVGmu(k)=sum(mu.*H)/H_total;
end
%figure(1)
%plot(CC,AVGmu-D)
%xlabel('K')
%ylabel('death')
%title('play')
%grid on

% [SRvalues',Hx',CC',VHR',num_strains']%%%%%%%%%%%%%%%

figure(2)
plot(Hx,Vx)
xlabel('H TOTAL')
ylabel('V TOTAL')
grid on

figure(3)
plot(Hx,AVGmu)
xlabel('H TOTAL')
ylabel('AVERAGE INTRINSIC GROWTH RATE')
grid on

figure(4)
plot(Hx,ShannonDiversity)
xlabel('H TOTAL')
ylabel('SHANNON DIVERSITY')
grid on

%figure(5)
%plot(Hx,AVGmu-D)
%xlabel('H TOTAL')
%ylabel('DEATH RATE BY PHAGE')
%grid on

%figure(6)
%plot(AVGmu, VHR)
%xlabel('Average Intrinsic Growth Rate')
%ylabel('VMR')
%grid on

%figure(9)
%plot(Hx./HT, AVGmu-D)
%xlabel('N/K')
%ylabel('Death')
%grid on



% **************
function s=sigp(sigma,n)
o=ones(n);
sig=sigma.^[0:n-1];
s=zeros(n);
for k=0:n-1
    v=diag(o,k)*sig(k+1);
    s=s+diag(v,k);
end
function r=rhop(rho,n)
rh=rho.^[0:n-1];
r=repmat(rh,n,1);
r=triu(r);
function b=beta(rho,sigma,n)
r=rhop(rho,n);
s=sigp(sigma,n);
b=r.*s;



