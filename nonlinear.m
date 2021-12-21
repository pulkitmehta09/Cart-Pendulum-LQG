function dX = nonlinear(t,X,F)
m1 = 100;
m2 = 100;
M = 1000;
l1 = 20;
l2 = 10;
g = 9.81;
x = X(1);
dx = X(2);
th1 = X(3);
dth1 = X(4);
th2 = X(5);
dth2 = X(6);
dX = zeros(6,1);
dX(1) = dx;
dX(2) = (F-m1*(g*sin(th1)*cos(th1)+l1*sin(th1)*(dth1)^2)-m2*(g*sin(th2)*cos(th2)+l2*sin(th2)*(dth2)^2))/(M+m1*(sin(th1))^2+m2*(sin(th2))^2);
dX(3) = dth1;
dX(4) = ((cos(th1))/l1)*(F-m1*(g*sin(th1)*cos(th1)+l1*sin(th1)*(dth1)^2)-m2*sin(th2)*cos(th2)+l2*sin(th2)*(dth2)^2)/(M+m1*(sin(th1))^2+m2*(sin(th2))^2)-(g*sin(th1))/l1;
dX(5) = dth2;
dX(6) = ((cos(th2))/l2)*(F-m1*(g*sin(th1)*cos(th1)+l1*sin(th1)*(dth1)^2)-m2*sin(th2)*cos(th2)+l2*sin(th2)*(dth2)^2)/(M+m1*(sin(th1))^2+m2*(sin(th2))^2)-(g*sin(th2))/l2;

end
