% Code RK4:
clf
t_end=150;
h=1/200;
t=0:h:t_end;
% step size
n=floor(t_end/h);
P2=h/2;

%% Parameters values
% p     = 0.03;    % Immunity loss
% Lambda= 15;      % Birth/recruitment rate (not too low)
% beta  = 0.0005;  % Infection rate (LOWERED!)
% d     = 0.01;    % Natural death rate
% gama  = 0.7;     % Recovery coefficient
% mu    = 0.6;     % Recovery rate
% b     = 0.01;    % Saturation incidence
% delta = 0.01;    % Disease-induced death rate
% m     = 0.01;    % Migration/other removal


p     = 0.03;    % Immunity loss
Lambda= 1;      % Birth/recruitment rate (not too low)
beta  = 0.001;  % Infection rate (LOWERED!)
d     = 0.01;    % Natural death rate
gama  = 0.4;     % Recovery coefficient
mu    = 0.6;     % Recovery rate
b     = 0.01;    % Saturation incidence
delta = 0.01;    % Disease-induced death rate
m     = 0.01;    % Migration/other removal


%% Classes
S0= 990; I0= 10;  R0=0;
%%%%
S10= 1100; I10= 20; R10=5;
%%%
S20= 1200; I20= 30;  R20=10;
%%%%
S30= 1300; I30= 40;  R30=15;
%%%
S40= 1400; I40= 50; R40=20;

%%%
S= zeros(1,n+1);  I= zeros(1,n+1);   R= zeros(1,n+1);
%%%
S1= zeros(1,n+1); I1= zeros(1,n+1);  R1= zeros(1,n+1);

S2= zeros(1,n+1); I2= zeros(1,n+1);  R2= zeros(1,n+1);

S3= zeros(1,n+1); I3= zeros(1,n+1);  R3= zeros(1,n+1);

S4= zeros(1,n+1); I4= zeros(1,n+1);  R4= zeros(1,n+1);
%%%%%%%%%%
S(1)=S0; I(1)=I0;  R(1)=R0;

S1(1)=S10; I1(1)=I10;  R1(1)=R10;

S2(1)=S20; I2(1)=I20; R2(1)=R20;

S3(1)=S30; I3(1)=I30; R3(1)=R30;

S4(1)=S40; I4(1)=I40; R4(1)=R40;
%%

for i =1:n  
     k1_1 =Lambda-beta*S(i)*I(i)-d*S(i)+p*R(i);
     k2_1 =beta*S(i)*I(i)-(gama*mu*I(i)/(1+b*mu*I(i)))-d*I(i)-delta*I(i)-m*I(i);
     k3_1 = (gama*mu*I(i)/(1+b*mu*I(i)))-d*R(i)-p*R(i)+m*I(i);
     
     k1_2 =Lambda-beta*(S(i)+P2*k1_1)*(I(i)+P2*k2_1)-d*(S(i)+P2*k1_1)+p*(R(i)+P2*k3_1);
     k2_2 =beta*(S(i)+P2*k1_1)*(I(i)+P2*k2_1)-(gama*mu*(I(i)+P2*k2_1)/(1+b*mu*(I(i)+P2*k2_1)))-(d+delta+m)*(I(i)+P2*k2_1);
     k3_2 = (gama*mu*(I(i)+P2*k2_1)/(1+b*mu*(I(i)+P2*k2_1)))-(d+p)*(R(i)+P2*k3_1)+m*(I(i)+P2*k2_1);
     
     k1_3 =Lambda-beta*(S(i)+P2*k1_2)*(I(i)+P2*k2_2)-d*(S(i)+P2*k1_2)+p*(R(i)+P2*k3_2);
     k2_3 =beta*(S(i)+P2*k1_2)*(I(i)+P2*k2_2)-(gama*mu*(I(i)+P2*k2_2)/(1+b*mu*(I(i)+P2*k2_2)))-(d+delta+m)*(I(i)+P2*k2_2);
     k3_3 = (gama*mu*(I(i)+P2*k2_2)/(1+b*mu*(I(i)+P2*k2_2)))-(d+p)*(R(i)+P2*k3_2)+m*(I(i)+P2*k2_2);
     
     k1_4 =Lambda-beta*(S(i)+P2*k1_3)*(I(i)+P2*k2_3)-d*(S(i)+P2*k1_3)+p*(R(i)+P2*k3_3);
     k2_4 =beta*(S(i)+P2*k1_3)*(I(i)+P2*k2_3)-(gama*mu*(I(i)+P2*k2_3)/(1+b*mu*(I(i)+P2*k2_3)))-(d+delta+m)*(I(i)+P2*k2_3);
     k3_4 = (gama*mu*(I(i)+P2*k2_3)/(1+b*mu*(I(i)+P2*k2_3)))-(d+p)*(R(i)+P2*k3_3)+m*(I(i)+P2*k2_3);
    
    %  S
    S(i+1) = S(i)+(h/6)*(k1_1+2*k1_2+2*k1_3+k1_4);
    % Infected I  
    I(i+1) = I(i)+(h/6)*(k2_1+2*k2_2+2*k2_3+k2_4);
    % R
    R(i+1) = R(i)+(h/6)*(k3_1+2*k3_2+2*k3_3+k3_4);
    end

    for i =1:n  
     k1_1 =Lambda-beta*S1(i)*I1(i)-d*S1(i)+p*R1(i);
     k2_1 =beta*S1(i)*I1(i)-(gama*mu*I1(i)/(1+b*mu*I1(i)))-d*I1(i)-delta*I1(i)-m*I1(i);
     k3_1 = (gama*mu*I1(i)/(1+b*mu*I1(i)))-d*R1(i)-p*R1(i)+m*I1(i);
     
     k1_2 =Lambda-beta*(S1(i)+P2*k1_1)*(I1(i)+P2*k2_1)-d*(S1(i)+P2*k1_1)+p*(R1(i)+P2*k3_1);
     k2_2 =beta*(S1(i)+P2*k1_1)*(I1(i)+P2*k2_1)-(gama*mu*(I1(i)+P2*k2_1)/(1+b*mu*(I1(i)+P2*k2_1)))-(d+delta+m)*(I1(i)+P2*k2_1);
     k3_2 = (gama*mu*(I1(i)+P2*k2_1)/(1+b*mu*(I1(i)+P2*k2_1)))-(d+p)*(R1(i)+P2*k3_1)+m*(I1(i)+P2*k2_1);
     
     k1_3 =Lambda-beta*(S1(i)+P2*k1_2)*(I1(i)+P2*k2_2)-d*(S1(i)+P2*k1_2)+p*(R1(i)+P2*k3_2);
     k2_3 =beta*(S1(i)+P2*k1_2)*(I1(i)+P2*k2_2)-(gama*mu*(I1(i)+P2*k2_2)/(1+b*mu*(I1(i)+P2*k2_2)))-(d+delta+m)*(I1(i)+P2*k2_2);
     k3_3 = (gama*mu*(I1(i)+P2*k2_2)/(1+b*mu*(I1(i)+P2*k2_2)))-(d+p)*(R1(i)+P2*k3_2)+m*(I1(i)+P2*k2_2);
     
     k1_4 =Lambda-beta*(S1(i)+P2*k1_3)*(I1(i)+P2*k2_3)-d*(S1(i)+P2*k1_3)+p*(R1(i)+P2*k3_3);
     k2_4 =beta*(S1(i)+P2*k1_3)*(I1(i)+P2*k2_3)-(gama*mu*(I1(i)+P2*k2_3)/(1+b*mu*(I1(i)+P2*k2_3)))-(d+delta+m)*(I1(i)+P2*k2_3);
     k3_4 = (gama*mu*(I1(i)+P2*k2_3)/(1+b*mu*(I1(i)+P2*k2_3)))-(d+p)*(R1(i)+P2*k3_3)+m*(I1(i)+P2*k2_3);
    
    %  S
    S1(i+1) = S1(i)+(h/6)*(k1_1+2*k1_2+2*k1_3+k1_4);
    % Infected I  
    I1(i+1) = I1(i)+(h/6)*(k2_1+2*k2_2+2*k2_3+k2_4);
    % R
    R1(i+1) = R1(i)+(h/6)*(k3_1+2*k3_2+2*k3_3+k3_4);
    end
    
    for i =1:n  
     k1_1 =Lambda-beta*S2(i)*I2(i)-d*S2(i)+p*R2(i);
     k2_1 =beta*S2(i)*I2(i)-(gama*mu*I2(i)/(1+b*mu*I2(i)))-d*I2(i)-delta*I2(i)-m*I2(i);
     k3_1 = (gama*mu*I2(i)/(1+b*mu*I2(i)))-d*R2(i)-p*R2(i)+m*I2(i);
     
     k1_2 =Lambda-beta*(S2(i)+P2*k1_1)*(I2(i)+P2*k2_1)-d*(S2(i)+P2*k1_1)+p*(R2(i)+P2*k3_1);
     k2_2 =beta*(S2(i)+P2*k1_1)*(I2(i)+P2*k2_1)-(gama*mu*(I2(i)+P2*k2_1)/(1+b*mu*(I2(i)+P2*k2_1)))-(d+delta+m)*(I2(i)+P2*k2_1);
     k3_2 = (gama*mu*(I2(i)+P2*k2_1)/(1+b*mu*(I2(i)+P2*k2_1)))-(d+p)*(R2(i)+P2*k3_1)+m*(I2(i)+P2*k2_1);
     
     k1_3 =Lambda-beta*(S2(i)+P2*k1_2)*(I2(i)+P2*k2_2)-d*(S2(i)+P2*k1_2)+p*(R2(i)+P2*k3_2);
     k2_3 =beta*(S2(i)+P2*k1_2)*(I2(i)+P2*k2_2)-(gama*mu*(I2(i)+P2*k2_2)/(1+b*mu*(I2(i)+P2*k2_2)))-(d+delta+m)*(I2(i)+P2*k2_2);
     k3_3 = (gama*mu*(I2(i)+P2*k2_2)/(1+b*mu*(I2(i)+P2*k2_2)))-(d+p)*(R2(i)+P2*k3_2)+m*(I2(i)+P2*k2_2);
     
     k1_4 =Lambda-beta*(S2(i)+P2*k1_3)*(I2(i)+P2*k2_3)-d*(S2(i)+P2*k1_3)+p*(R2(i)+P2*k3_3);
     k2_4 =beta*(S2(i)+P2*k1_3)*(I2(i)+P2*k2_3)-(gama*mu*(I2(i)+P2*k2_3)/(1+b*mu*(I2(i)+P2*k2_3)))-(d+delta+m)*(I2(i)+P2*k2_3);
     k3_4 = (gama*mu*(I2(i)+P2*k2_3)/(1+b*mu*(I2(i)+P2*k2_3)))-(d+p)*(R2(i)+P2*k3_3)+m*(I2(i)+P2*k2_3);
    
    %  S
    S2(i+1) = S2(i)+(h/6)*(k1_1+2*k1_2+2*k1_3+k1_4);
    % Infected I  
    I2(i+1) = I2(i)+(h/6)*(k2_1+2*k2_2+2*k2_3+k2_4);
    % R
    R2(i+1) = R2(i)+(h/6)*(k3_1+2*k3_2+2*k3_3+k3_4);
    end
    
    
     for i =1:n  
     k1_1 =Lambda-beta*S3(i)*I3(i)-d*S3(i)+p*R3(i);
     k2_1 =beta*S3(i)*I3(i)-(gama*mu*I3(i)/(1+b*mu*I3(i)))-d*I3(i)-delta*I3(i)-m*I3(i);
     k3_1 = (gama*mu*I3(i)/(1+b*mu*I3(i)))-d*R3(i)-p*R3(i)+m*I3(i);
     
     k1_2 =Lambda-beta*(S3(i)+P2*k1_1)*(I3(i)+P2*k2_1)-d*(S3(i)+P2*k1_1)+p*(R3(i)+P2*k3_1);
     k2_2 =beta*(S3(i)+P2*k1_1)*(I3(i)+P2*k2_1)-(gama*mu*(I3(i)+P2*k2_1)/(1+b*mu*(I3(i)+P2*k2_1)))-(d+delta+m)*(I3(i)+P2*k2_1);
     k3_2 = (gama*mu*(I3(i)+P2*k2_1)/(1+b*mu*(I3(i)+P2*k2_1)))-(d+p)*(R3(i)+P2*k3_1)+m*(I3(i)+P2*k2_1);
     
     k1_3 =Lambda-beta*(S3(i)+P2*k1_2)*(I3(i)+P2*k2_2)-d*(S3(i)+P2*k1_2)+p*(R3(i)+P2*k3_2);
     k2_3 =beta*(S3(i)+P2*k1_2)*(I3(i)+P2*k2_2)-(gama*mu*(I3(i)+P2*k2_2)/(1+b*mu*(I3(i)+P2*k2_2)))-(d+delta+m)*(I3(i)+P2*k2_2);
     k3_3 = (gama*mu*(I3(i)+P2*k2_2)/(1+b*mu*(I3(i)+P2*k2_2)))-(d+p)*(R3(i)+P2*k3_2)+m*(I3(i)+P2*k2_2);
     
     k1_4 =Lambda-beta*(S3(i)+P2*k1_3)*(I3(i)+P2*k2_3)-d*(S3(i)+P2*k1_3)+p*(R3(i)+P2*k3_3);
     k2_4 =beta*(S3(i)+P2*k1_3)*(I3(i)+P2*k2_3)-(gama*mu*(I3(i)+P2*k2_3)/(1+b*mu*(I3(i)+P2*k2_3)))-(d+delta+m)*(I3(i)+P2*k2_3);
     k3_4 = (gama*mu*(I3(i)+P2*k2_3)/(1+b*mu*(I3(i)+P2*k2_3)))-(d+p)*(R3(i)+P2*k3_3)+m*(I3(i)+P2*k2_3);
    
    %  S
    S3(i+1) = S3(i)+(h/6)*(k1_1+2*k1_2+2*k1_3+k1_4);
    % Infected I  
    I3(i+1) = I3(i)+(h/6)*(k2_1+2*k2_2+2*k2_3+k2_4);
    % R
    R3(i+1) = R3(i)+(h/6)*(k3_1+2*k3_2+2*k3_3+k3_4);
     end
    
 for i =1:n  
     
     k1_1 =Lambda-beta*S4(i)*I4(i)-d*S4(i)+p*R4(i);
     k2_1 =beta*S4(i)*I4(i)-(gama*mu*I4(i)/(1+b*mu*I4(i)))-d*I4(i)-delta*I4(i)-m*I4(i);
     k3_1 = (gama*mu*I4(i)/(1+b*mu*I4(i)))-d*R4(i)-p*R4(i)+m*I4(i);
     
     k1_2 =Lambda-beta*(S4(i)+P2*k1_1)*(I4(i)+P2*k2_1)-d*(S4(i)+P2*k1_1)+p*(R4(i)+P2*k3_1);
     k2_2 =beta*(S4(i)+P2*k1_1)*(I4(i)+P2*k2_1)-(gama*mu*(I4(i)+P2*k2_1)/(1+b*mu*(I4(i)+P2*k2_1)))-(d+delta+m)*(I4(i)+P2*k2_1);
     k3_2 = (gama*mu*(I4(i)+P2*k2_1)/(1+b*mu*(I4(i)+P2*k2_1)))-(d+p)*(R4(i)+P2*k3_1)+m*(I4(i)+P2*k2_1);
     
     k1_3 =Lambda-beta*(S4(i)+P2*k1_2)*(I4(i)+P2*k2_2)-d*(S4(i)+P2*k1_2)+p*(R4(i)+P2*k3_2);
     k2_3 =beta*(S4(i)+P2*k1_2)*(I4(i)+P2*k2_2)-(gama*mu*(I4(i)+P2*k2_2)/(1+b*mu*(I4(i)+P2*k2_2)))-(d+delta+m)*(I4(i)+P2*k2_2);
     k3_3 = (gama*mu*(I4(i)+P2*k2_2)/(1+b*mu*(I4(i)+P2*k2_2)))-(d+p)*(R4(i)+P2*k3_2)+m*(I4(i)+P2*k2_2);
     
     k1_4 =Lambda-beta*(S4(i)+P2*k1_3)*(I4(i)+P2*k2_3)-d*(S4(i)+P2*k1_3)+p*(R4(i)+P2*k3_3);
     k2_4 =beta*(S4(i)+P2*k1_3)*(I4(i)+P2*k2_3)-(gama*mu*(I4(i)+P2*k2_3)/(1+b*mu*(I4(i)+P2*k2_3)))-(d+delta+m)*(I4(i)+P2*k2_3);
     k3_4 = (gama*mu*(I4(i)+P2*k2_3)/(1+b*mu*(I4(i)+P2*k2_3)))-(d+p)*(R4(i)+P2*k3_3)+m*(I4(i)+P2*k2_3);
    
    %  S
    S4(i+1) = S4(i)+(h/6)*(k1_1+2*k1_2+2*k1_3+k1_4);
    % Infected I  
    I4(i+1) = I4(i)+(h/6)*(k2_1+2*k2_2+2*k2_3+k2_4);
    % R
    R4(i+1) = R4(i)+(h/6)*(k3_1+2*k3_2+2*k3_3+k3_4);
    end
%% Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
plot(t,S,'b:',t,S1,'r-.',t,S2,'k--',t,S3,'g',t,S4,'m:','LineWidth',2.2);
legend('Case1','Case2','Case3','Case4','Case5')
xlabel('time t (Days)')
ylabel('Susceptible class')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
plot(t,I,'b:',t,I1,'r-.',t,I2,'k--',t,I3,'g',t,I4,'m:','LineWidth',2.2);
legend('Case1','Case2','Case3','Case4','Case5')
xlabel('time t (Days)')
ylabel('Infected class')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3)
plot(t,R,'b:',t,R1,'r-.',t,R2,'k--',t,R3,'g',t,R4,'m:','LineWidth',2.2);
legend('Case1','Case2','Case3','Case4','Case5')
xlabel('time t (Days)')
ylabel('Recoevered class')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
save Alldata