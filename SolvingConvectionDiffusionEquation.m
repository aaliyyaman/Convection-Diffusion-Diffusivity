
clear
close all
clc


L=10000; % in ft
dT=1; %in day
dX=50;
v=1.5; % velocity of conc. in ft/day
D=93.01; % diffusion coef. in ft^2/day
Cinitial=0;
Cb=1; % boundary concentration

alpha=D*dT/dX^2;
a=alpha;

beta=v*dT/dX;
b=beta;


% EXPLİCİT
A=diag([1-2*b-3*a,(1-b-2*a)*ones(1,18),1-b-a])+diag(((a+b)*ones(1,19)),-1)+diag((a*ones(1,19)),1);
cn=zeros(20,1)*Cinitial;
Q=zeros(20,1);
Q(1)=-2*Cb*(b+a);

C_exp=zeros(500,20);
for i=1:500
    Cnp1=A*cn+Q;
    cn=Cnp1;
    C_exp(i,:)=abs(Cnp1');
end


% IMPLICIT
AA=diag([1+2*b+3*a,(1+b+2*a)*ones(1,18),1+b+a])+diag(((-b-a)*ones(1,19)),-1)+diag((-a*ones(1,19)),1);
cn=zeros(20,1);
Q=zeros(20,1);
Q(1)=2*Cb*(a+b);

C_imp=zeros(500,20);
for i=1:500
    B=cn+Q;
    Cnp1=AA\B;
    cn=Cnp1;
    C_imp(i,:)=Cnp1';
end

%ANALYTİC SOLUTİON 
cc=linspace(0,1000,21);
t=200;
AnSol=1/2*(erfc((cc-v*t)/(2*sqrt(D*t)))+exp(v*cc/D).*erfc((cc+v*t)/(2*sqrt(D*t))));



% PLOTTİNG RESULTS
cex=[1 C_exp(200,:)];
cim=[1 C_imp(200,:)];
figure
hold on
plot(cc, cex, '+', 'MarkerSize', 8);

plot(cc, cim, 'o', 'MarkerSize', 6);

plot(cc, AnSol, 'k*', 'MarkerSize', 6);

xlabel('Delta X, from x=0 to x=1000 ft')
ylabel('Convection, C')


legend('Explicit', 'Implicit','Analytical Solution')

title('Comparing solutions for Convection Diffusiton')

hold off

