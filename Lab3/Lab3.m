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
eps = 0.000001;
Id = inf;
numIter = 0;
yGained = zeros(M,N);
