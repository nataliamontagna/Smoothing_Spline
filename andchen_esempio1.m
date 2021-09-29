n=150;
a=-0.2;
b=0.4;

sigma=0.2;
x=a+sort(rand(n,1))*(b-a);
f=@(t) cos(2*pi*t)+0.3*sin(10*pi*t)+0.2*t;
y=f(x)+sigma*randn(n,1);

p=2;
gamma=(b-a)^(2*p-1);
[Ut,Vt]=generators1((x-a)./(b-a),p);

% lambda=1.e-10
lambda1=1.e-10;

subplot(1,3,1)
plot(x,y,'o')
hold on
subplot(1,3,1)
plot(x,f(x),'--k')
hold on

[Wt,z]=potrf1(Ut,Vt,n*lambda1/gamma);
[ctilde,d]=smoothing_spline_reg(Ut,Wt,z,y);

[phi,zeta]=andersen_chen(a,b,p,n,x);
fhat=@(t) 0;
for k=1:p
    fhat=@(t) fhat(t)+d(k)*phi{k}(t);
end
for j=1:n
    fhat=@(t) fhat(t)+(ctilde(j))/gamma * zeta{j}(t);
end

t=linspace(a,b,1000);
subplot(1,3,1)
plot(t,fhat(t),'-r')
title('\lambda = 10^{-10}');

% lambda=1.e-7
lambda2=1.e-7;
subplot(1,3,2)
plot(x,y,'o')
hold on
subplot(1,3,2)
plot(x,f(x),'--k')
hold on

[Wt,z]=potrf1(Ut,Vt,n*lambda2/gamma);
[ctilde,d]=smoothing_spline_reg(Ut,Wt,z,y);
[phi,zeta]=andersen_chen(a,b,p,n,x);
fhat=@(t) 0;
for k=1:p
    fhat=@(t) fhat(t)+d(k)*phi{k}(t);
end
for j=1:n
    fhat=@(t) fhat(t)+(ctilde(j))/gamma * zeta{j}(t);
end
t=linspace(a,b,1000);
subplot(1,3,2)
plot(t,fhat(t),'-r')
title('\lambda = 10^{-7}');

% lambda=1.e-3
lambda3=1.e-3;

subplot(1,3,3)
plot(x,y,'o')
hold on
subplot(1,3,3)
plot(x,f(x),'--k')
hold on

[Wt,z]=potrf1(Ut,Vt,n*lambda3/gamma);
[ctilde,d]=smoothing_spline_reg(Ut,Wt,z,y);

[phi,zeta]=andersen_chen(a,b,p,n,x);
fhat=@(t) 0;
for k=1:p
    fhat=@(t) fhat(t)+d(k)*phi{k}(t);
end
for j=1:n
    fhat=@(t) fhat(t)+(ctilde(j))/gamma * zeta{j}(t);
end

t=linspace(a,b,1000);
subplot(1,3,3)
plot(t,fhat(t),'-r')
title('\lambda = 10^{-3}');