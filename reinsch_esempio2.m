Y=load('mcyc.dat');
tobs=Y(:,1); yobs=Y(:,2);
[t,y,w]=splav(tobs,yobs);

N=length(t);
n=N-1;
a=0;
b=60;
deltay = sqrt(N) * ones(N,1);
D = diag(deltay);

% S=100
S1 = 100;
norm(D\y)^2-S1

subplot(1,3,1)
plot(tobs,yobs,'o')
hold on

[a1,b1,c1,d1] = reinsch(t,y,S1,deltay);

for i = 1:n
    x = t(i):0.1:t(i+1);
    subplot(1,3,1)
    plot(x,a1(i)+b1(i)*(x-t(i))+c1(i)*(x-t(i)).^2+d1(i)*(x-t(i)).^3,'-r');
    hold on
end
title('S=100');

% S=300
S2=300;
norm(D\y)^2-S2

subplot(1,3,2)
plot(t,y,'o')
hold on

[a1,b1,c1,d1] = reinsch(t,y,S2,deltay);

for i = 1:n
    x = t(i):0.1:t(i+1);
    subplot(1,3,2)
    plot(x,a1(i)+b1(i)*(x-t(i))+c1(i)*(x-t(i)).^2+d1(i)*(x-t(i)).^3,'-r');
    hold on
end
title('S=300');

% S=600
S3=600;
norm(D\y)^2-S3

subplot(1,3,3)
plot(t,y,'o')
hold on

[a1,b1,c1,d1] = reinsch(t,y,S3,deltay);

for i = 1:n
    x = t(i):0.1:t(i+1);
    subplot(1,3,3)
    plot(x,a1(i)+b1(i)*(x-t(i))+c1(i)*(x-t(i)).^2+d1(i)*(x-t(i)).^3,'-r');
    hold on
end
title('S=600');