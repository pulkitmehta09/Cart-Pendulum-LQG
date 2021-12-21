syms x x_d t1 t1_d t2 t2_d t F M m1 m2 g l1 l2
f1 = x_d;
f2 = (F-m1*(g*sin(t1)*cos(t1)+l1*sin(t1)*(t1_d)^2)-m2*(g*sin(t2)*cos(t2)+l2*sin(t2)*(t2_d)^2))/(M+m1*(sin(t1))^2+m2*(sin(t2))^2);
f3 = t1_d;
f4 = ((cos(t1))/l1)*(F-m1*(g*sin(t1)*cos(t1)+l1*sin(t1)*(t1_d)^2)-m2*sin(t2)*cos(t2)+l2*sin(t2)*(t2_d)^2)/(M+m1*(sin(t1))^2+m2*(sin(t2))^2)-(g*sin(t1))/l1;
f5 = t2_d; 
f6 = ((cos(t2))/l2)*(F-m1*(g*sin(t1)*cos(t1)+l1*sin(t1)*(t1_d)^2)-m2*sin(t2)*cos(t2)+l2*sin(t2)*(t2_d)^2)/(M+m1*(sin(t1))^2+m2*(sin(t2))^2)-(g*sin(t2))/l2;

J = jacobian([f1,f2,f3,f4,f5,f6],[x,x_d,t1,t1_d,t2,t2_d]);
% Matrix A
A = subs(J,[x,x_d,t1,t1_d,t2,t2_d],[0,0,0,0,0,0])
J_u = jacobian([f1,f2,f3,f4,f5,f6],[F]);
B = subs(J_u,[x,x_d,t1,t1_d,t2,t2_d],[0,0,0,0,0,0])