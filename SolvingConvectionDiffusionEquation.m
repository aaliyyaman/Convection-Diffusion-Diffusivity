
clear
close all
clc



dT=1; %in day
dX=50;
v=1.5; % velocity of conc. in ft/day
D=93.01; % diffusion coef. in ft^2/day
Cinitial=0;
Cb=1; % boundary concentration

L=1000;
a=D*dT/dX^2;
b=v*dT/dX;

ngrid=25;
time=[500 250 100 50];


for tt=1:size(time,2)
    aa=time(tt);
    calculation(aa,a,b,L,v,D,Cinitial,Cb,ngrid);
end


function calculation(time,a,b,L,v,D,Cinitial,Cb,ngrid)
    
    
    
    % EXPLİCİT
    A=diag([1-2*b-3*a,(1-b-2*a)*ones(1,ngrid-2),1-b-a])+diag(((a+b)*ones(1,ngrid-1)),-1)+diag((a*ones(1,ngrid-1)),1);
    cn=zeros(ngrid,1)*Cinitial;
    Q=zeros(ngrid,1);
    Q(1)=-2*Cb*(b+a);
    
    C_exp=zeros(time,ngrid);
    for i=1:time
        Cnp1=A*cn+Q;
        cn=Cnp1;
        C_exp(i,:)=abs(Cnp1');
    end
    
    
    % IMPLICIT
    AA=diag([1+2*b+3*a,(1+b+2*a)*ones(1,ngrid-2),1+b+a])+diag(((-b-a)*ones(1,ngrid-1)),-1)+diag((-a*ones(1,ngrid-1)),1);
    cn=zeros(ngrid,1)*Cinitial;
    Q=zeros(ngrid,1);
    Q(1)=2*Cb*(a+b);
    
    C_imp=zeros(time,ngrid);
    for i=1:time
        B=cn+Q;
        Cnp1=AA\B;
        cn=Cnp1;
        C_imp(i,:)=Cnp1';
    end
    
    %ANALYTİC SOLUTİON
    cc=linspace(0,L,ngrid+1);
    t=time;
    AnSol=1/2*(erfc((cc-v*t)/(2*sqrt(D*t)))+exp(v*cc/D).*erfc((cc+v*t)/(2*sqrt(D*t))));
    
    
    
    % PLOTTİNG RESULTS
    cex=[1 C_exp(time,:)];
    cim=[1 C_imp(time,:)];
    
    figure
    hold on
    plot(cc, cex, '+', 'MarkerSize', 8);
    
    plot(cc, cim, 'o', 'MarkerSize', 6);
    
    plot(cc, AnSol, 'k*', 'MarkerSize', 6);
    
    xlabel(sprintf('Delta X, from x=0 to x=%sft',num2str(L)))
    ylabel('Convection, C')
    
    
    legend('Explicit', 'Implicit','Analytical Solution')
    
    title(sprintf(' Ploting for %sth day ' , num2str(time)));
    
    hold off
    
end
