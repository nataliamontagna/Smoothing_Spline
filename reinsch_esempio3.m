load carbon12alpha

N=length(angle);
n=N-1;
a=0;
b=4.2;

x=angle;
y=counts;

deltay = sqrt(N) * ones(N,1);
D = diag(deltay);

% S=60
S1 = 60;
norm(D\y)^2-S1

subplot(1,3,1)
plot(x,y,'o')
hold on

[a1,b1,c1,d1] = reinsch(x,y,S1,deltay);

for i = 1:n
    t = x(i):0.001:x(i+1);
    subplot(1,3,1)
    plot(t,a1(i)+b1(i)*(t-x(i))+c1(i)*(t-x(i)).^2+d1(i)*(t-x(i)).^3,'-r');
    hold on
end
axis([a b -50 350]);
title('S=60');

% S=400
S2 = 400;
norm(D\y)^2-S2

subplot(1,3,2)
plot(x,y,'o')
hold on

[a1,b1,c1,d1] = reinsch(x,y,S2,deltay);

for i = 1:n
    t = x(i):0.001:x(i+1);
    subplot(1,3,2)
    plot(t,a1(i)+b1(i)*(t-x(i))+c1(i)*(t-x(i)).^2+d1(i)*(t-x(i)).^3,'-r');
    hold on
end
axis([a b -50 350]);
title('S=400');

% S=1000
S3 = 1000;
norm(D\y)^2-S3

subplot(1,3,3)
plot(x,y,'o')
hold on

[a1,b1,c1,d1] = reinsch(x,y,S3,deltay);

for i = 1:n
    t = x(i):0.001:x(i+1);
    subplot(1,3,3)
    plot(t,a1(i)+b1(i)*(t-x(i))+c1(i)*(t-x(i)).^2+d1(i)*(t-x(i)).^3,'-r');
    hold on
end
axis([a b -50 350]);
title('S=1000');