load carbon12alpha

n=length(angle);
a=0;
b=4.2;

x=angle;
y=counts;

p=2;
gamma=(b-a)^(2*p-1);
[Ut,Vt]=generators1((x-a)./(b-a),p);

% lambda=1.e-7
lambda1=1.e-7;

subplot(1,3,1)
plot(x,y,'o')
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
plot(t,fhat(t))
axis([a b -50 350]);
title('\lambda = 10^{-7}')

% lambda=1.e-5
lambda2=1.e-5;
subplot(1,3,2)
plot(x,y,'o')
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

subplot(1,3,2)
plot(t,fhat(t))
axis([a b -50 350]);
title('\lambda = 10^{-5}');

% lambda=1.e-3
lambda3=1.e-3;

subplot(1,3,3)
plot(x,y,'o')
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

subplot(1,3,3)
plot(t,fhat(t))
axis([a b -50 350]);
title('\lambda = 10^{-3}')