Y=load('mcyc.dat');
tobs=Y(:,1); yobs=Y(:,2);
[t,y,w]=splav(tobs,yobs);

N=length(t);
n=N-1;
a=0;
b=60;

deltay = sqrt(N) * ones(N,1);

S = 300;

plot(tobs,yobs,'o')
hold on

[a1,b1,c1,d1,q] = reinsch(t,y,S,deltay);

s=[];
sval=[];

for i = 1:n
    x = t(i):0.1:t(i+1);
    s=[s, x];
    xval=a1(i)+b1(i)*(x-t(i))+c1(i)*(x-t(i)).^2+d1(i)*(x-t(i)).^3;
    plot(x,xval,'-r');
    hold on
    sval=[sval, xval];
end

lambda = 1/(2*q);

p=2;
gamma=(b-a)^(2*p-1);
[Ut,Vt]=generators1((t-a)./(b-a),p);
[Wt,z]=potrf1(Ut,Vt,N*lambda/gamma);
[ctilde,d]=smoothing_spline_reg(Ut,Wt,z,y);
[phi,zeta]=andersen_chen(a,b,p,N,t);
fhat=@(t) 0;
for k=1:p
    fhat=@(t) fhat(t)+d(k)*phi{k}(t);
end
for j=1:N
    fhat=@(t) fhat(t)+(ctilde(j))/gamma * zeta{j}(t);
end
sval_ac=fhat(s);
plot(s,sval_ac,'--b')

err=norm(sval-sval_ac)/norm(sval)
