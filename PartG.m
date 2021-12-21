clear all, close all, clc

% Variable definitions

m1 = 100;
m2 = 100;
M = 1000;
l1 = 20;
l2 = 10;
g = 9.81;
dt = 0.1;
t = dt:dt:100;

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
D = zeros(size(C1,1),size(B,2));
sys1 = ss(A,B,C1,D);

%% LQR Controller
Q = [1 0 0 0 0 0; 0 10 0 0 0 0; 0 0 2000 0 0 0; 0 0 0 2000 0 0; 0 0 0 0 2000 0; 0 0 0 0 0 2000]; 
R = 0.001;
K = lqr(A,B,Q,R);
Ac = [(A-B*K)];
Bc = [B];
Cc = [C1];
Dc = [D];

sys_cl = ss(Ac,Bc,Cc,Dc);

%% Luenberger Observer

Bd = 0.1*eye(6);
Vn = 0.1;
l1 = (lqr(A',C1',Bd,Vn))';
sysL1 = ss(A-(l1*C1),[B l1],C1,0);

%% LQG applied to Non-linear system

[t,out] = ode45(@(t,x)nonlinearObserver(t,x,C1,-K*x,l1),t,X_0);
figure()
hold on
plot(t,out(:,1))
ylabel('State variable')
xlabel('time (sec)')
title('LQG response to Non-linear system for output vector: x(t)')
legend('x')
hold off
