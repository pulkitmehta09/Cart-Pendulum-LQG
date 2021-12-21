clear all, close all, clc

% Variable definitions
m1 = 100;
m2 = 100;
M = 1000;
l1 = 20;
l2 = 10;
g = 9.81;
 

%% Observability Check
A = [0 1 0 0 0 0; 0 0 -m1*g/M 0 -m2*g/M 0; 0 0 0 1 0 0; 0 0 -((M*g)+(m1*g))/(M*l1) 0 -g*m2/(M*l1) 0; 0 0 0 0 0 1; 0 0 -m1*g/(M*l2) 0 -((M*g)+(m2*g))/(M*l2) 0];
B = [0; 1/M; 0; 1/(l1*M); 0; 1/(l2*M)];
D = [0;0;0];
% Case 1: Output vector: x
C1 = [1 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];
Obs1 = [C1' (C1*A)' (C1*A^2)' (C1*A^3)' (C1*A^4)' (C1*A^5)'];
rank_Obs1 = rank(Obs1);
fprintf('Rank of controllability matrix with output vector as x: %d\n',rank_Obs1)

% Case 2 Output vectors: theta1, theta2
C2 = [0 0 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 1 0];
Obs2 = [C2' (C2*A)' (C2*A^2)' (C2*A^3)' (C2*A^4)' (C2*A^5)'];
rank_Obs2 = rank(Obs2);
fprintf('Rank of controllability matrix with output vectors as theta1 and theta2: %d\n',rank_Obs2)

% Case 3 Output vectors: x, theta2
C3 = [1 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 1 0];
Obs3 = [C3' (C3*A)' (C3*A^2)' (C3*A^3)' (C3*A^4)' (C3*A^5)'];
rank_Obs3 = rank(Obs3);
fprintf('Rank of controllability matrix with output vector as x and theta2: %d\n',rank_Obs3)


% Case 2 Output vectors: x, theta1, theta2
C4 = [1 0 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 1 0];
Obs4 = [C4' (C4*A)' (C4*A^2)' (C4*A^3)' (C4*A^4)' (C4*A^5)'];
rank_Obs4 = rank(Obs4);
fprintf('Rank of controllability matrix with output vector as x, theta1 and theta2: %d\n',rank_Obs4)

