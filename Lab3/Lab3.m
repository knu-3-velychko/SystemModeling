clear
clc
yExperimental = dlmread('y9.txt',' ');
plot(0:0.2:50, yExperimental(1,:));
[M,N] = size(yExperimental);
c(1) = 0.14;
c(2) = 0.1;
c(3) = 0.2;
c(4) = 0.12;
m(1) = 11;
m(2) = 28;
m(3) = 23;
h = 0.2;
t = 0.0;
tk = 50.0;
eps = 0.000001;
I = inf;
numIter = 0;
yGained = zeros(M,N);

function res = getA(m, c)
    A = zeros(6,6);
    A(1, 2) = 1;
    A(2, 1) = -(c(1)+c(2))/m(1);
    A(2, 3) = c(2)/m(1);
    A(3, 4) = 1;
    A(4, 1) = c(2)/m(2);
    A(4, 3) = -(c(2)+c(3))/m(2);
    A(4, 5) = c(3)/m(2);
    A(5, 6) = 1;
    A(6, 3) = c(3)/m(3);
    A(6, 5) = -(c(3)+c(4))/m(3);
end

function res = fy(A, y)
  res = A*y
end

function res = fU(A, U, y)
  C = zeros(6,3);
  C(2, 1) = (y(3)-y(1))/m(1);
  C(2, 2) = -(c(2)*y(3) - (c(1)+c(2))*y(3))/(m(1)*m(1));
  C(4, 1) = (y(1) - y(3))/m(2);
  C(6, 3) = -(c(3)*y(3) - (c(3) + c(4))*y(5))/(m(3)*m(3));
  res = A*U + C
end

function res = RungeKuttU(U, h, y, m, c)
  k1 = h * fU(U, y, m, c);
  k2 = h * fU(U + k1/2.0, y, m, c);
  k3 = h * fU(U + k2/2.0, y, m, c);
  k4 = h * fU(U + k3, y, m, c);
  res = U + (k1 + 2*k2 + 2*k3 + k4)/6.0;
endfunction

function res = RungeKuttY(y, h, m, c)
  k1 = h * fy(y, m, c);
  k2 = h * fy(y + k1/2.0, m, c);
  k3 = h * fy(y + k2/2.0, m, c);
  k4 = h * fy(y + k3, m, c);
  res = y + (k1 + 2*k2 + 2*k3 + k4)/6.0;
endfunction


while (I > eps)
  numIter = numIter+1;
  y = yExperimental(:,1);
  
end