close all;
clear;
%initial guess
x0=[1000 -1000 120 120];
options=optimset('Algorithm','active-set','TolFun',1e-9);
[x,fval]=fmincon(@Obj_22_2a,x0,[],[],[],[],[],[],@Con_22_2a,options)

function f = Obj_22_2a(x)
    %Compute geometry
    L1=(x(1)^2+1000^2)^0.5;
    L2=(x(2)^2+1000^2)^0.5;
    A1=pi*x(3)^2/4;
    A2=pi*x(4)^2/4;
    rho=7800;
    %Compute mass
    f=rho*(A1*L1+A2*L2)/1e9;
end

function [c, ceq]=Con_22_2a(x)
    %Compute geometry, angles are relative to global x-axis
    theta1=atan2(x(1),-1000);
    theta2=atan2(x(2),-1000);
    L1=(x(1)^2+1000^2)^0.5;
    L2=(x(2)^2+1000^2)^0.5;
    A1=pi*x(3)^2/4;
    A2=pi*x(4)^2/4;
    I1=pi*x(3)^4/64;
    I2=pi*x(4)^4/64;
    %Compute forces, 1st eq: sumFy=0, 2nd eq: sumFx=0
    M=[sin(theta1) sin(theta2); cos(theta1) cos(theta2)];
    Y=[10000 0]';
    F=inv(M)*Y;
    %Compute stress
    sigma1=abs(F(1))/A1;
    sigma2=abs(F(2))/A2;
    %Compute buckling forces
    Fbck1=pi^2*2.1e5*I1/L1^2;
    Fbck2=pi^2*2.1e5*I2/L2^2;
    %Inequality constraints
    c(1)=2-x(3);
    c(2)=2-x(4);
    c(3)=sigma1-230/1.5;
    c(4)=sigma2-230/1.5;
    %If tension force, then compare 0.0 with buckling force
    %If compression force, then compare absolute value with buckling force
    Fcmp1=min(0,F(1));
    Fcmp2=min(0,F(2));
    c(5)=abs(Fcmp1)-Fbck1/2.5;
    c(6)=abs(Fcmp2)-Fbck2/2.5;
    ceq=[];
end