load carbon12alpha

n=length(angle);
a=0;
b=4.2;

x=angle;
y=counts;

p=2;
gamma=(b-a)^(2*p-1);
[Ut,Vt]=generators1((x-a)./(b-a),p);

% GML
q1=fminbnd(@(v) l_gml(Ut,Vt,n,gamma,y,v), -10, 10);
lambda1=10^q1;
subplot(1,2,1)
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
subplot(1,2,1)
plot(t,fhat(t))
axis([a b -50 350]);
title(['\lambda = ' num2str(lambda1)])

% GCV
q2=fminbnd(@(v) l_gcv(Ut,Vt,n,gamma,y,v), -10, 10);
lambda2=10^q2;
subplot(1,2,2)
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
subplot(1,2,2)
plot(t,fhat(t))
axis([a b -50 350]);
title(['\lambda = ' num2str(lambda2)])