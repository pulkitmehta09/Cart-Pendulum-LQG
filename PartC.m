clear all, close all, clc

%% Controllability
syms m1 g m2 M l1 l2

A = [0 1 0 0 0 0; 0 0 -m1*g/M 0 -m2*g/M 0; 0 0 0 1 0 0; 0 0 -((M*g)+(m1*g))/(M*l1) 0 -g*m2/(M*l1) 0; 0 0 0 0 0 1; 0 0 -m1*g/(M*l2) 0 -((M*g)+(m2*g))/(M*l2) 0];
B = [0; 1/M; 0; 1/(l1*M); 0; 1/(l2*M)];
C = [1 0 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 1 0];
D = [0;0;0];

Co = [B A*B (A^2)*B (A^3)*B (A^4)*B (A^5)*B];
rank = rank(Co);
fprintf('Rank of controllability matrix: %d\n',rank)

det_Co = det(Co);
disp('Determinant of controllability matrix:')
disp(det_Co)
