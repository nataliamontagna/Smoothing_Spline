n = 150;
a = -0.2;
b = 0.4;

sigma = 0.2;
x = a + sort(rand(n,1)) * (b-a);
f = @(t) cos(2*pi*t)+0.3*sin(10*pi*t)+0.2*t;
y = f(x) + sigma * randn(n,1);
deltay = sigma * ones(n,1);
D=diag(deltay);

% S=70
S1=70;
norm(D\y)^2-S1

subplot(1,3,1)
plot(x,y,'o')
hold on
subplot(1,3,1)
plot(x,f(x),'--k')
hold on

[a1,b1,c1,d1] = reinsch(x,y,S1,deltay);

for i = 1:n-1
    t = x(i):0.0001:x(i+1);
    subplot(1,3,1)
    plot(t,a1(i)+b1(i)*(t-x(i))+c1(i)*(t-x(i)).^2+d1(i)*(t-x(i)).^3,'-r');
    hold on
end
title('S=70');

% S=150
S2 = 150;
norm(D\y)^2-S2

subplot(1,3,2)
plot(x,y,'o')
hold on
subplot(1,3,2)
plot(x,f(x),'--k')
hold on

[a1,b1,c1,d1] = reinsch(x,y,S2,deltay);

for i = 1:n-1
    t = x(i):0.0001:x(i+1);
    subplot(1,3,2)
    plot(t,a1(i)+b1(i)*(t-x(i))+c1(i)*(t-x(i)).^2+d1(i)*(t-x(i)).^3,'-r');
    hold on
end
title('S=150');

% S=270
S3 = 270;
norm(D\y)^2-S3

subplot(1,3,3)
plot(x,y,'o')
hold on
subplot(1,3,3)
plot(x,f(x),'--k')
hold on

[a1,b1,c1,d1] = reinsch(x,y,S3,deltay);

for i = 1:n-1
    t = x(i):0.0001:x(i+1);
    subplot(1,3,3)
    plot(t,a1(i)+b1(i)*(t-x(i))+c1(i)*(t-x(i)).^2+d1(i)*(t-x(i)).^3,'-r');
    hold on
end
title('S=270');