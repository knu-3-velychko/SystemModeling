clear
clc
y = dlmread('f19.txt');
dt = 0.01;
T = 5;
t = 0:dt:T;

plot(t,y),grid;
f = fft(y);

n = length(y);

c1(1:n) = 0;
c2(1:n) = 0;
c(1:n) = 0;

for k = 1:n
    cos_sum = 0;
    sin_sum = 0;
    for  m = 1:n
        cos_sum = cos_sum + y(m)*cos(2*pi*(k-1)*(m-1)/n);
        sin_sum = sin_sum + y(m)*sin(2*pi*(k-1)*(m-1)/n); 
    end
    c1(k) = cos_sum/n;
    c2(k) = sin_sum/n;
    c(k) = sqrt(c1(k)*c1(k)+c2(k)*c2(k));
end

figure;
plot(t,c),grid;

figure;
plot(abs(f)),grid;

for k=3: n/2
    
   if isequal(c(k), max(c(k-2:k+2)))
       disp(k);
       index=k; 
       break;
   end   
end

index=index-1
f=index/T

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
disp(A);
   
C=[sum(y.*t.^3);
   sum(y.*t.^2);
   sum(y.*t);
   sum(y.*sin(t.*2*pi*f));
   sum(y)];
disp(C);


x=A\C;

approximation = x(1)*t.^3+x(2)*t.^2+x(3).*t+x(4)*sin(2*pi*f*t)+x(5);
figure;
plot(t,approximation),grid;
err = 0;
for i = 1:n 
    err = err + (approximation(i)-y(i))*(approximation(i)-y(i));
end;
err