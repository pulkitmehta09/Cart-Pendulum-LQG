clear all, close all, clc

% Variable definitions
m1 = 100;
m2 = 100;
M = 1000;
l1 = 20;
l2 = 10;
g = 9.81;
dt = 0.1;
tspan = dt:dt:100;
% Initial conditions
x_0 = 1;
theta1_0 = deg2rad(10);
theta2_0 = deg2rad(20);
X_0 = [1 0 theta1_0 0 theta2_0 0];

%% Observability Check
A = [0 1 0 0 0 0; 0 0 -m1*g/M 0 -m2*g/M 0; 0 0 0 1 0 0; 0 0 -((M*g)+(m1*g))/(M*l1) 0 -g*m2/(M*l1) 0; 0 0 0 0 0 1; 0 0 -m1*g/(M*l2) 0 -((M*g)+(m2*g))/(M*l2) 0];
B = [0; 1/M; 0; 1/(l1*M); 0; 1/(l2*M)];

% Case 1: Output vector: x
C1 = [1 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];
Obs1 = obsv(A,C1);
rank_Obs1 = rank(Obs1);
fprintf('Rank of controllability matrix with output vector as x: %d\n',rank_Obs1)

% Case 2 Output vectors: theta1, theta2
C2 = [0 0 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 1 0];
Obs2 = obsv(A,C2);
rank_Obs2 = rank(Obs2);
fprintf('Rank of controllability matrix with output vectors as theta1 and theta2: %d\n',rank_Obs2)

% Case 3 Output vectors: x, theta2
C3 = [1 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 1 0];
Obs3 = obsv(A,C3);
rank_Obs3 = rank(Obs3);
fprintf('Rank of controllability matrix with output vector as x and theta2: %d\n',rank_Obs3)

% Case 2 Output vectors: x, theta1, theta2
C4 = [1 0 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 1 0];
Obs4 = obsv(A,C4);
rank_Obs4 = rank(Obs4);
fprintf('Rank of controllability matrix with output vector as x, theta1 and theta2: %d\n',rank_Obs4)

D = zeros(size(C1,1),size(B,2));
sys1 = ss(A,B,C1,D);
sys3 = ss(A,B,C3,D);
sys4 = ss(A,B,C4,D);

%% Luenberger Observer

Bd = 0.1*eye(6);
Vn = 0.1;

% [l1,P,E] = lqe(A,Bd,C1,Bd,Vn*eye(3));
l1 = (lqr(A',C1',Bd,Vn))';
sysL1 = ss(A-(l1*C1),[B l1],C1,0);
l3 = (lqr(A',C3',Bd,Vn))';
sysL3 = ss(A-(l3*C3),[B l3],C3,0);
l4 = (lqr(A',C4',Bd,Vn))';
sysL4 = ss(A-(l4*C4),[B l4],C4,0);


u = 0*tspan;
u(100:length(tspan)) = 1;

[y1,t] = lsim(sys1,u,tspan);
[x1,t] = lsim(sysL1,[u;y1'],tspan);

[y3,t] = lsim(sys3,u,tspan);
[x3,t] = lsim(sysL3,[u;y3'],tspan);

[y4,t] = lsim(sys4,u,tspan);
[x4,t] = lsim(sysL4,[u;y4'],tspan);


figure(1)
hold on
plot(t,y1(:,1),'g',LineWidth=2);
plot(t,x1(:,1),'k--','LineWidth',1.0)
ylabel('State Variables')
xlabel('Time(sec)')
legend('x(t)','Estimated x(t)')
title('Response for output vector at step input: x(t)')
hold off

figure(2)
hold on
plot(t,y3(:,1),'g',LineWidth=2);
plot(t,x3(:,1),'k--','LineWidth',1.0)
plot(t,y3(:,3),'r',LineWidth=2);
plot(t,x3(:,3),'k--','LineWidth',1.0)
ylabel('State Variables')
xlabel('Time(sec)')
legend('x(t)','theta2(t)','Estimated x(t)','Estimated theta2(t)')
title('Response for output vector at step input: x(t),theta2(t)')
hold off

figure(3)
hold on
plot(t,y3(:,1),'g',LineWidth=2);
plot(t,x3(:,1),'k--','LineWidth',1.0)
plot(t,y3(:,2),'c',LineWidth=2);
plot(t,x3(:,2),'k--','LineWidth',1.0)
plot(t,y3(:,3),'r',LineWidth=2);
plot(t,x3(:,3),'k--','LineWidth',1.0)
ylabel('State Variables')
xlabel('Time(sec)')
legend('x(t)','theta1(t)','theta2(t)','Estimated x(t)','Estimated theta1(t)','Estimated theta2(t)')
title('Response for output vector at step input: x(t),theta1(t),theta2(t)')
hold off

%% Linear Model response to Observer

[t,l_out1] = ode45(@(t,x)linearObserver(t,x,C1,l1),tspan,X_0);
figure(4)
hold on
plot(t,l_out1(:,1))
ylabel('state variables')
xlabel('time (sec)')
title('Linear system Observer for output vector: x(t)')
legend('x')
hold off

[t,l_out3] = ode45(@(t,x)linearObserver(t,x,C3,l3),tspan,X_0);
figure(5)
hold on
plot(t,l_out3(:,1))
plot(t,l_out3(:,5))
ylabel('state variables')
xlabel('time (sec)')
title('Linear system Observer for output vectors: x(t),theta2(t)')
legend('x','theta2')
hold off
 
[t,l_out4] = ode45(@(t,x)linearObserver(t,x,C4,l4),tspan,X_0);
figure(6)
hold on
plot(t,l_out4(:,1))
plot(t,l_out4(:,3))
plot(t,l_out4(:,5))
ylabel('state variables')
xlabel('time (sec)')
title('Linear system Observer for output vectors: x(t),theta1(t),theta2(t)')
legend('x','theta1','theta2')
hold off
 
%% Non-Linear Model response to Observer

[t,nl_out1] = ode45(@(t,x)nonlinearObserver(t,x,C1,1,l1),tspan,X_0);
figure(7)
hold on
plot(t,nl_out1(:,1))
ylabel('state variables')
xlabel('time (sec)')
title('Non Linear system Observer for output vector: x(t)')
legend('x')
hold off

[t,nl_out3] = ode45(@(t,x)nonlinearObserver(t,x,C3,1,l3),tspan,X_0);
figure(8)
hold on
plot(t,nl_out3(:,1))
plot(t,nl_out3(:,5))
ylabel('state variables')
xlabel('time (sec)')
title('Non-Linear system Observer for output vectors: x(t),theta2(t)')
legend('x','theta2')
hold off
 
[t,nl_out4] = ode45(@(t,x)nonlinearObserver(t,x,C4,1,l4),tspan,X_0);
figure(9)
hold on
plot(t,nl_out4(:,1))
plot(t,nl_out4(:,3))
plot(t,nl_out4(:,5))
ylabel('state variables')
xlabel('time (sec)')
title('Non-Linear system Observer for output vectors: x(t),theta1(t),theta2(t)')
legend('x','theta1','theta2')
hold off


