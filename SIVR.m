% Code RK4:
clf
t_end=150;
h=1/200;
t=0:h:t_end;
% step size
n=floor(t_end/h);
P2=h/2;

%% Parameters values
% Lambda=0.9; beta=0.001; p=0.5; alpha=0.3; d=0.09; mu=0.005; q=0.6; gama=0.3;
% b=0.09; delta=0.5; m=0.04; kapa=0.3;

% p= 0.4; Lambda = 10; beta = 0.009; alpha = 0.1; q=0.02; d= 0.01; gama  = 0.1; mu= 0.2; b= 0.05;
% delta  = 0.01; m= 0.01; kapa= 0.05;

p= 0.4; Lambda =0.9; beta = 0.009; alpha = 0.1; q=0.09; d= 0.05; gama  = 0.3; mu= 0.2; b= 0.05;
delta  = 0.01; m= 0.01; kapa= 0.05;

%% Classes
S0= 990; I0= 10; V0=0; R0=0;
%%%%
S10= 1100; I10= 20; V10=10; R10=5;
%%%
S20= 1200; I20= 30; V20=20; R20=10;
%%%%
S30= 1300; I30= 40; V30=30; R30=15;
%%%
S40= 1400; I40= 50; V40=40; R40=20;

%%%
S= zeros(1,n+1);  I= zeros(1,n+1);  V= zeros(1,n+1); R= zeros(1,n+1);
%%%
S1= zeros(1,n+1); I1= zeros(1,n+1); V1= zeros(1,n+1); R1= zeros(1,n+1);

S2= zeros(1,n+1); I2= zeros(1,n+1); V2= zeros(1,n+1); R2= zeros(1,n+1);

S3= zeros(1,n+1); I3= zeros(1,n+1); V3= zeros(1,n+1); R3= zeros(1,n+1);

S4= zeros(1,n+1); I4= zeros(1,n+1); V4= zeros(1,n+1); R4= zeros(1,n+1);
%%%%%%%%%%
S(1)=S0; I(1)=I0; V(1)=V0; R(1)=R0;

S1(1)=S10; I1(1)=I10; V1(1)=V10; R1(1)=R10;

S2(1)=S20; I2(1)=I20; V2(1)=V20; R2(1)=R20;

S3(1)=S30; I3(1)=I30; V3(1)=V30; R3(1)=R30;

S4(1)=S40; I4(1)=I40; V4(1)=V40; R4(1)=R40;
%%
for i =1:n  
     k1_1 =(1-p)*Lambda-(beta*S(i)*I(i)/(1+alpha*I(i)))-d*S(i)+q*V(i);
     k2_1 =(beta*S(i)*I(i)/(1+alpha*I(i)))-(gama*mu*I(i)/(1+b*mu*I(i)))-d*I(i)-delta*I(i)-m*I(i);
     k3_1 =p*Lambda+m*I(i)-d*V(i)-q*V(i)-kapa*V(i); 
     k4_1 = kapa*V(i)+(gama*mu*I(i)/(1+b*mu*I(i)))-d*R(i);
     
     k1_2 =(1-p)*Lambda-(beta*(S(i)+P2*k1_1)*(I(i)+P2*k2_1)/(1+alpha*(I(i)+P2*k2_1)))-d*(S(i)+P2*k1_1)+q*(V(i)+P2*k3_1);
     k2_2 =(beta*(S(i)+P2*k1_1)*(I(i)+P2*k2_1)/(1+alpha*(I(i)+P2*k2_1)))-(gama*mu*(I(i)+P2*k2_1)/(1+b*mu*(I(i)+P2*k2_1)))-(d+delta+m)*(I(i)+P2*k2_1);
     k3_2 =p*Lambda+m*(I(i)+P2*k2_1)-(d+q+kapa)*(V(i)+P2*k3_1); 
     k4_2 = kapa*(V(i)+P2*k3_1)+(gama*mu*(I(i)+P2*k2_1)/(1+b*mu*(I(i)+P2*k2_1)))-d*(R(i)+P2*k4_1);
    
     k1_3 =(1-p)*Lambda-(beta*(S(i)+P2*k1_2)*(I(i)+P2*k2_2)/(1+alpha*(I(i)+P2*k2_2)))-d*(S(i)+P2*k1_2)+q*(V(i)+P2*k3_2);
     k2_3 =(beta*(S(i)+P2*k1_2)*(I(i)+P2*k2_2)/(1+alpha*(I(i)+P2*k2_2)))-(gama*mu*(I(i)+P2*k2_2)/(1+b*mu*(I(i)+P2*k2_2)))-(d+delta+m)*(I(i)+P2*k2_2);
     k3_3 =p*Lambda+m*(I(i)+P2*k2_2)-(d+q+kapa)*(V(i)+P2*k3_2); 
     k4_3 = kapa*(V(i)+P2*k3_2)+(gama*mu*(I(i)+P2*k2_2)/(1+b*mu*(I(i)+P2*k2_2)))-d*(R(i)+P2*k4_2);
     
     k1_4 =(1-p)*Lambda-(beta*(S(i)+P2*k1_3)*(I(i)+P2*k2_3)/(1+alpha*(I(i)+P2*k2_3)))-d*(S(i)+P2*k1_3)+q*(V(i)+P2*k3_3);
     k2_4 =(beta*(S(i)+P2*k1_3)*(I(i)+P2*k2_3)/(1+alpha*(I(i)+P2*k2_3)))-(gama*mu*(I(i)+P2*k2_3)/(1+b*mu*(I(i)+P2*k2_3)))-(d+delta+m)*(I(i)+P2*k2_3);
     k3_4 =p*Lambda+m*(I(i)+P2*k2_3)-(d+q+kapa)*(V(i)+P2*k3_3); 
     k4_4 = kapa*(V(i)+P2*k3_3)+(gama*mu*(I(i)+P2*k2_3)/(1+b*mu*(I(i)+P2*k2_3)))-d*(R(i)+P2*k4_3);
    
    %  S
    S(i+1) = S(i)+(h/6)*(k1_1+2*k1_2+2*k1_3+k1_4);
    % Infected I  
    I(i+1) = I(i)+(h/6)*(k2_1+2*k2_2+2*k2_3+k2_4);
    % V
    V(i+1) = V(i)+(h/6)*(k3_1+2*k3_2+2*k3_3+k3_4);
    % R
    R(i+1) = R(i)+(h/6)*(k4_1+2*k4_2+2*k4_3+k4_4);
    end

    for i =1:n  
     k1_1 =(1-p)*Lambda-(beta*S1(i)*I1(i)/(1+alpha*I1(i)))-d*S1(i)+q*V1(i);
     k2_1 =(beta*S1(i)*I1(i)/(1+alpha*I1(i)))-(gama*mu*I1(i)/(1+b*mu*I1(i)))-(d+delta+m)*I1(i);
     k3_1 =p*Lambda+m*I1(i)-(d+q+kapa)*V1(i); 
     k4_1 = kapa*V1(i)+(gama*mu*I1(i)/(1+b*mu*I1(i)))-d*R1(i);
     
     k1_2 =(1-p)*Lambda-(beta*(S1(i)+P2*k1_1)*(I1(i)+P2*k2_1)/(1+alpha*(I1(i)+P2*k2_1)))-d*(S1(i)+P2*k1_1)+q*(V1(i)+P2*k3_1);
     k2_2 =(beta*(S1(i)+P2*k1_1)*(I1(i)+P2*k2_1)/(1+alpha*(I1(i)+P2*k2_1)))-(gama*mu*(I1(i)+P2*k2_1)/(1+b*mu*(I1(i)+P2*k2_1)))-(d+delta+m)*(I1(i)+P2*k2_1);
     k3_2 =p*Lambda+m*(I1(i)+P2*k2_1)-(d+q+kapa)*(V1(i)+P2*k3_1); 
     k4_2 = kapa*(V1(i)+P2*k3_1)+(gama*mu*(I1(i)+P2*k2_1)/(1+b*mu*(I1(i)+P2*k2_1)))-d*(R1(i)+P2*k4_1);
    
     k1_3 =(1-p)*Lambda-(beta*(S1(i)+P2*k1_2)*(I1(i)+P2*k2_2)/(1+alpha*(I1(i)+P2*k2_2)))-d*(S1(i)+P2*k1_2)+q*(V1(i)+P2*k3_2);
     k2_3 =(beta*(S1(i)+P2*k1_2)*(I1(i)+P2*k2_2)/(1+alpha*(I1(i)+P2*k2_2)))-(gama*mu*(I1(i)+P2*k2_2)/(1+b*mu*(I1(i)+P2*k2_2)))-(d+delta+m)*(I1(i)+P2*k2_2);
     k3_3 =p*Lambda+m*(I1(i)+P2*k2_2)-(d+q+kapa)*(V1(i)+P2*k3_2); 
     k4_3 = kapa*(V1(i)+P2*k3_2)+(gama*mu*(I1(i)+P2*k2_2)/(1+b*mu*(I1(i)+P2*k2_2)))-d*(R1(i)+P2*k4_2);
     
     k1_4 =(1-p)*Lambda-(beta*(S1(i)+P2*k1_3)*(I1(i)+P2*k2_3)/(1+alpha*(I1(i)+P2*k2_3)))-d*(S1(i)+P2*k1_3)+q*(V1(i)+P2*k3_3);
     k2_4 =(beta*(S1(i)+P2*k1_3)*(I1(i)+P2*k2_3)/(1+alpha*(I1(i)+P2*k2_3)))-(gama*mu*(I1(i)+P2*k2_3)/(1+b*mu*(I1(i)+P2*k2_3)))-(d+delta+m)*(I1(i)+P2*k2_3);
     k3_4 =p*Lambda+m*(I1(i)+P2*k2_3)-(d+q+kapa)*(V1(i)+P2*k3_3); 
     k4_4 = kapa*(V1(i)+P2*k3_3)+(gama*mu*(I1(i)+P2*k2_3)/(1+b*mu*(I1(i)+P2*k2_3)))-d*(R1(i)+P2*k4_3);
    
    %  S
    S1(i+1) = S1(i)+(h/6)*(k1_1+2*k1_2+2*k1_3+k1_4);
    % Infected I  
    I1(i+1) = I1(i)+(h/6)*(k2_1+2*k2_2+2*k2_3+k2_4);
    % V
    V1(i+1) = V1(i)+(h/6)*(k3_1+2*k3_2+2*k3_3+k3_4);
    % R
    R1(i+1) = R1(i)+(h/6)*(k4_1+2*k4_2+2*k4_3+k4_4);
    end
    
    for i =1:n  
     k1_1 =(1-p)*Lambda-(beta*S2(i)*I2(i)/(1+alpha*I2(i)))-d*S2(i)+q*V2(i);
     k2_1 =(beta*S2(i)*I2(i)/(1+alpha*I2(i)))-(gama*mu*I2(i)/(1+b*mu*I2(i)))-(d+delta+m)*I2(i);
     k3_1 =p*Lambda+m*I2(i)-(d+q+kapa)*V2(i); 
     k4_1 = kapa*V2(i)+(gama*mu*I2(i)/(1+b*mu*I2(i)))-d*R2(i);
     
     k1_2 =(1-p)*Lambda-(beta*(S2(i)+P2*k1_1)*(I2(i)+P2*k2_1)/(1+alpha*(I2(i)+P2*k2_1)))-d*(S2(i)+P2*k1_1)+q*(V2(i)+P2*k3_1);
     k2_2 =(beta*(S2(i)+P2*k1_1)*(I2(i)+P2*k2_1)/(1+alpha*(I2(i)+P2*k2_1)))-(gama*mu*(I2(i)+P2*k2_1)/(1+b*mu*(I2(i)+P2*k2_1)))-(d+delta+m)*(I2(i)+P2*k2_1);
     k3_2 =p*Lambda+m*(I2(i)+P2*k2_1)-(d+q+kapa)*(V2(i)+P2*k3_1); 
     k4_2 = kapa*(V2(i)+P2*k3_1)+(gama*mu*(I2(i)+P2*k2_1)/(1+b*mu*(I2(i)+P2*k2_1)))-d*(R2(i)+P2*k4_1);
    
     k1_3 =(1-p)*Lambda-(beta*(S2(i)+P2*k1_2)*(I2(i)+P2*k2_2)/(1+alpha*(I2(i)+P2*k2_2)))-d*(S2(i)+P2*k1_2)+q*(V2(i)+P2*k3_2);
     k2_3 =(beta*(S2(i)+P2*k1_2)*(I2(i)+P2*k2_2)/(1+alpha*(I2(i)+P2*k2_2)))-(gama*mu*(I2(i)+P2*k2_2)/(1+b*mu*(I2(i)+P2*k2_2)))-(d+delta+m)*(I2(i)+P2*k2_2);
     k3_3 =p*Lambda+m*(I2(i)+P2*k2_2)-(d+q+kapa)*(V2(i)+P2*k3_2); 
     k4_3 = kapa*(V2(i)+P2*k3_2)+(gama*mu*(I2(i)+P2*k2_2)/(1+b*mu*(I2(i)+P2*k2_2)))-d*(R2(i)+P2*k4_2);
     
     k1_4 =(1-p)*Lambda-(beta*(S2(i)+P2*k1_3)*(I2(i)+P2*k2_3)/(1+alpha*(I2(i)+P2*k2_3)))-d*(S2(i)+P2*k1_3)+q*(V2(i)+P2*k3_3);
     k2_4 =(beta*(S2(i)+P2*k1_3)*(I2(i)+P2*k2_3)/(1+alpha*(I2(i)+P2*k2_3)))-(gama*mu*(I2(i)+P2*k2_3)/(1+b*mu*(I2(i)+P2*k2_3)))-(d+delta+m)*(I2(i)+P2*k2_3);
     k3_4 =p*Lambda+m*(I2(i)+P2*k2_3)-(d+q+kapa)*(V2(i)+P2*k3_3); 
     k4_4 = kapa*(V2(i)+P2*k3_3)+(gama*mu*(I2(i)+P2*k2_3)/(1+b*mu*(I2(i)+P2*k2_3)))-d*(R2(i)+P2*k4_3);
    
    %  S
    S2(i+1) = S2(i)+(h/6)*(k1_1+2*k1_2+2*k1_3+k1_4);
    % Infected I  
    I2(i+1) = I2(i)+(h/6)*(k2_1+2*k2_2+2*k2_3+k2_4);
    % V
    V2(i+1) = V2(i)+(h/6)*(k3_1+2*k3_2+2*k3_3+k3_4);
    % R
    R2(i+1) = R2(i)+(h/6)*(k4_1+2*k4_2+2*k4_3+k4_4);
    end
    
    
     for i =1:n  
     k1_1 =(1-p)*Lambda-(beta*S3(i)*I3(i)/(1+alpha*I3(i)))-d*S3(i)+q*V3(i);
     k2_1 =(beta*S3(i)*I3(i)/(1+alpha*I3(i)))-(gama*mu*I3(i)/(1+b*mu*I3(i)))-(d+delta+m)*I3(i);
     k3_1 =p*Lambda+m*I3(i)-(d+q+kapa)*V3(i); 
     k4_1 = kapa*V3(i)+(gama*mu*I3(i)/(1+b*mu*I3(i)))-d*R3(i);
     
     k1_2 =(1-p)*Lambda-(beta*(S3(i)+P2*k1_1)*(I3(i)+P2*k2_1)/(1+alpha*(I3(i)+P2*k2_1)))-d*(S3(i)+P2*k1_1)+q*(V3(i)+P2*k3_1);
     k2_2 =(beta*(S3(i)+P2*k1_1)*(I3(i)+P2*k2_1)/(1+alpha*(I3(i)+P2*k2_1)))-(gama*mu*(I3(i)+P2*k2_1)/(1+b*mu*(I3(i)+P2*k2_1)))-(d+delta+m)*(I3(i)+P2*k2_1);
     k3_2 =p*Lambda+m*(I3(i)+P2*k2_1)-(d+q+kapa)*(V3(i)+P2*k3_1); 
     k4_2 = kapa*(V3(i)+P2*k3_1)+(gama*mu*(I3(i)+P2*k2_1)/(1+b*mu*(I3(i)+P2*k2_1)))-d*(R3(i)+P2*k4_1);
    
     k1_3 =(1-p)*Lambda-(beta*(S3(i)+P2*k1_2)*(I3(i)+P2*k2_2)/(1+alpha*(I3(i)+P2*k2_2)))-d*(S3(i)+P2*k1_2)+q*(V3(i)+P2*k3_2);
     k2_3 =(beta*(S3(i)+P2*k1_2)*(I3(i)+P2*k2_2)/(1+alpha*(I3(i)+P2*k2_2)))-(gama*mu*(I3(i)+P2*k2_2)/(1+b*mu*(I3(i)+P2*k2_2)))-(d+delta+m)*(I3(i)+P2*k2_2);
     k3_3 =p*Lambda+m*(I3(i)+P2*k2_2)-(d+q+kapa)*(V3(i)+P2*k3_2); 
     k4_3 = kapa*(V3(i)+P2*k3_2)+(gama*mu*(I3(i)+P2*k2_2)/(1+b*mu*(I3(i)+P2*k2_2)))-d*(R3(i)+P2*k4_2);
     
     k1_4 =(1-p)*Lambda-(beta*(S3(i)+P2*k1_3)*(I3(i)+P2*k2_3)/(1+alpha*(I3(i)+P2*k2_3)))-d*(S3(i)+P2*k1_3)+q*(V3(i)+P2*k3_3);
     k2_4 =(beta*(S3(i)+P2*k1_3)*(I3(i)+P2*k2_3)/(1+alpha*(I3(i)+P2*k2_3)))-(gama*mu*(I3(i)+P2*k2_3)/(1+b*mu*(I3(i)+P2*k2_3)))-(d+delta+m)*(I3(i)+P2*k2_3);
     k3_4 =p*Lambda+m*(I3(i)+P2*k2_3)-(d+q+kapa)*(V3(i)+P2*k3_3); 
     k4_4 = kapa*(V3(i)+P2*k3_3)+(gama*mu*(I3(i)+P2*k2_3)/(1+b*mu*(I3(i)+P2*k2_3)))-d*(R3(i)+P2*k4_3);
    
    %  S
    S3(i+1) = S3(i)+(h/6)*(k1_1+2*k1_2+2*k1_3+k1_4);
    % Infected I  
    I3(i+1) = I3(i)+(h/6)*(k2_1+2*k2_2+2*k2_3+k2_4);
    % V
    V3(i+1) = V3(i)+(h/6)*(k3_1+2*k3_2+2*k3_3+k3_4);
    % R
    R3(i+1) = R3(i)+(h/6)*(k4_1+2*k4_2+2*k4_3+k4_4);
     end
    
 for i =1:n  
     k1_1 =(1-p)*Lambda-(beta*S4(i)*I4(i)/(1+alpha*I4(i)))-d*S4(i)+q*V4(i);
     k2_1 =(beta*S4(i)*I4(i)/(1+alpha*I4(i)))-(gama*mu*I4(i)/(1+b*mu*I4(i)))-(d+delta+m)*I4(i);
     k3_1 =p*Lambda+m*I4(i)-(d+q+kapa)*V4(i); 
     k4_1 = kapa*V4(i)+(gama*mu*I4(i)/(1+b*mu*I4(i)))-d*R4(i);
     
     k1_2 =(1-p)*Lambda-(beta*(S4(i)+P2*k1_1)*(I4(i)+P2*k2_1)/(1+alpha*(I4(i)+P2*k2_1)))-d*(S4(i)+P2*k1_1)+q*(V4(i)+P2*k3_1);
     k2_2 =(beta*(S4(i)+P2*k1_1)*(I4(i)+P2*k2_1)/(1+alpha*(I4(i)+P2*k2_1)))-(gama*mu*(I4(i)+P2*k2_1)/(1+b*mu*(I4(i)+P2*k2_1)))-(d+delta+m)*(I4(i)+P2*k2_1);
     k3_2 =p*Lambda+m*(I4(i)+P2*k2_1)-(d+q+kapa)*(V4(i)+P2*k3_1); 
     k4_2 = kapa*(V4(i)+P2*k3_1)+(gama*mu*(I4(i)+P2*k2_1)/(1+b*mu*(I4(i)+P2*k2_1)))-d*(R4(i)+P2*k4_1);
    
     k1_3 =(1-p)*Lambda-(beta*(S4(i)+P2*k1_2)*(I4(i)+P2*k2_2)/(1+alpha*(I4(i)+P2*k2_2)))-d*(S4(i)+P2*k1_2)+q*(V4(i)+P2*k3_2);
     k2_3 =(beta*(S4(i)+P2*k1_2)*(I4(i)+P2*k2_2)/(1+alpha*(I4(i)+P2*k2_2)))-(gama*mu*(I4(i)+P2*k2_2)/(1+b*mu*(I4(i)+P2*k2_2)))-(d+delta+m)*(I4(i)+P2*k2_2);
     k3_3 =p*Lambda+m*(I4(i)+P2*k2_2)-(d+q+kapa)*(V4(i)+P2*k3_2); 
     k4_3 = kapa*(V4(i)+P2*k3_2)+(gama*mu*(I4(i)+P2*k2_2)/(1+b*mu*(I4(i)+P2*k2_2)))-d*(R4(i)+P2*k4_2);
     
     k1_4 =(1-p)*Lambda-(beta*(S4(i)+P2*k1_3)*(I4(i)+P2*k2_3)/(1+alpha*(I4(i)+P2*k2_3)))-d*(S4(i)+P2*k1_3)+q*(V4(i)+P2*k3_3);
     k2_4 =(beta*(S4(i)+P2*k1_3)*(I4(i)+P2*k2_3)/(1+alpha*(I4(i)+P2*k2_3)))-(gama*mu*(I4(i)+P2*k2_3)/(1+b*mu*(I4(i)+P2*k2_3)))-(d+delta+m)*(I4(i)+P2*k2_3);
     k3_4 =p*Lambda+m*(I4(i)+P2*k2_3)-(d+q+kapa)*(V4(i)+P2*k3_3); 
     k4_4 = kapa*(V4(i)+P2*k3_3)+(gama*mu*(I4(i)+P2*k2_3)/(1+b*mu*(I4(i)+P2*k2_3)))-d*(R4(i)+P2*k4_3);
    
    %  S
    S4(i+1) = S4(i)+(h/6)*(k1_1+2*k1_2+2*k1_3+k1_4);
    % Infected I  
    I4(i+1) = I4(i)+(h/6)*(k2_1+2*k2_2+2*k2_3+k2_4);
    % V
    V4(i+1) = V4(i)+(h/6)*(k3_1+2*k3_2+2*k3_3+k3_4);
    % R
    R4(i+1) = R4(i)+(h/6)*(k4_1+2*k4_2+2*k4_3+k4_4);
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
plot(t,V,'b:',t,V1,'r-.',t,V2,'k--',t,V3,'g',t,V4,'m:','LineWidth',2.2);
legend('Case1','Case2','Case3','Case4','Case5')
xlabel('time t (Days)')
ylabel(' Vaccinated class')
%%%%%%%%%%%%%%%%%%%%
figure(4)
plot(t,R,'b:',t,R1,'r-.',t,R2,'k--',t,R3,'g',t,R4,'m:','LineWidth',2.2);
legend('Case1','Case2','Case3','Case4','Case5')
xlabel('time t (Days)')
ylabel('Recoevered class')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
save Alldata