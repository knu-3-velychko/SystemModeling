clear
clc
input = dlmread('f19.txt');
dt = 0.01;
T = 5;
t = 0:dt:T;

plot(t,input),grid;

f = fft(input);
n = length(input);

c1(1:n) = 0;
c2(1:n) = 0;
c(1:n) = 0;

for k = 1:n
    cosFrequency = 0;
    sinFrequncy = 0;
    for  m = 1:n
        cosFrequency = cosFrequency + input(m)*cos(2*pi*(k-1)*(m-1)/n);
        sinFrequncy = sinFrequncy + input(m)*sin(2*pi*(k-1)*(m-1)/n); 
    end
    c1(k) = cosFrequency/n;
    c2(k) = sinFrequncy/n;
    c(k) = sqrt(c1(k)*c1(k)+c2(k)*c2(k));
end

figure;
plot(t,c),grid;

figure;
plot(abs(f)),grid;

figure;
plot(t,c),grid;
hold on
plot(abs(f));

for k=3: n/2
   if isequal(c(k), max(c(k-2:k+2)))
       k
       c(k)
       index=k; 
       break;
   end   
end

index=index-1;
f=index/T;

f

TPow(1:6)=0;
for k=1:6
    TPow(k)=sum(t.^k);
end

SinPow(1:3)=0;
for k=1:3
    SinPow(k)=sum((t.^k).*sin(t.*2*pi*f));
end

A(1:5,1:5)=[TPow(6), TPow(5), TPow(4), SinPow(3), TPow(3);
            TPow(5), TPow(4) ,TPow(3), SinPow(2), TPow(2);
            TPow(4) ,TPow(3), TPow(2), SinPow(1), TPow(1);
            SinPow(3), SinPow(2), SinPow(1), sum(sin(2*pi*f*t).^2),sum(sin(2*pi*f*t)) ; 
            TPow(3) ,TPow(2), TPow(1),sum(sin(2*pi*f*t)), n];
A
   
C=[sum(input.*t.^3);
   sum(input.*t.^2);
   sum(input.*t);
   sum(input.*sin(t.*2*pi*f));
   sum(input)];
   
C

a=A\C;

a

approximation = a(1)*t.^3+a(2)*t.^2+a(3).*t+a(4)*sin(2*pi*f*t)+a(5);
figure;
plot(t,approximation),grid;
figure;
plot(t,approximation),grid;
hold on
plot(t,input);
err = 0;
for i = 1:n 
    err = err + (approximation(i)-input(i))*(approximation(i)-input(i));
end;
err