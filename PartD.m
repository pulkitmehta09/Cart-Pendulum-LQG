clear all, close all, clc

% Variable definitions
m1 = 100;
m2 = 100;
M = 1000;
l1 = 20;
l2 = 10;
g = 9.81;
timespan = 0:0.05:100;
% Initial conditions
x_0 = 1;
theta1_0 = deg2rad(10);
theta2_0 = deg2rad(20);
X_0 = [1 0 theta1_0 0 theta2_0 0];

%% Linearized model
A = [0 1 0 0 0 0; 0 0 -m1*g/M 0 -m2*g/M 0; 0 0 0 1 0 0; 0 0 -((M*g)+(m1*g))/(M*l1) 0 -g*m2/(M*l1) 0; 0 0 0 0 0 1; 0 0 -m1*g/(M*l2) 0 -((M*g)+(m2*g))/(M*l2) 0];
B = [0; 1/M; 0; 1/(l1*M); 0; 1/(l2*M)];
C = [1 0 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 1 0];
D = [0;0;0];
Co = ctrb(A,B);
rank = rank(Co);
fprintf('Rank of controllability matrix: %d\n',rank)

%% LQR Controller

Q = [1 0 0 0 0 0; 0 10 0 0 0 0; 0 0 2000 0 0 0; 0 0 0 2000 0 0; 0 0 0 0 2000 0; 0 0 0 0 0 2000]; 
R = 0.001;
K = lqr(A,B,Q,R);
Ac = [(A-B*K)];
Bc = [B];
Cc = [C];
Dc = [D];

states = {'x' 'x_dot' 'theta1' 'theta1_dot' 'theta2' 'theta2_dot'};
inputs = {'u'};
outputs = {'x'; 'theta1'; 'theta2'};

sys_cl = ss(Ac,Bc,Cc,Dc,'statename',states,'inputname',inputs,'outputname',outputs);
figure(1)
step(sys_cl,200);
eig_Ac = eig(Ac);
disp('Eigenvalues of A-BK:')
display(eig_Ac)


%% Linear Model response to LQR
[t,l_out] = ode45(@(t,x)linear(t,x,-K*x),timespan,X_0);
figure(2);
hold on
plot(t,l_out(:,1))
plot(t,l_out(:,3))
plot(t,l_out(:,5))
ylabel('state variables')
xlabel('time (sec)')
title('LQR controller applied to linear system')
legend('x','theta1','theta2')


%% Non Linear Model response to LQR
[t,nl_out] = ode45(@(t,x)nonlinear(t,x,-K*x),timespan,X_0);
figure(3);
hold on
plot(t,nl_out(:,1))
plot(t,nl_out(:,3))
plot(t,nl_out(:,5))
ylabel('state variables')
xlabel('time (sec)')
title('LQR controller applied to non-linear system')
legend('x','theta1','theta2')


