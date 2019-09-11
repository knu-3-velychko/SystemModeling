data=dlmread('f19.txt','');
dt=0.01
T=5
t=0:dt:T;
plot(t,data);
xlabel("T");
ylabel("y(ti)");
title("Input data");