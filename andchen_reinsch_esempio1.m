n = 150;
a = -0.2;
b = 0.4;

sigma = 0.2;
x = a + sort(rand(n,1))*(b-a);
f = @(t) cos(2*pi*t) + 0.3*sin(10*pi*t) + 0.2*t;
y = f(x) + sigma * randn(n,1);
deltay = sqrt(n) * ones(n,1);

S=0.04;

plot(x,y,'o')
hold on
plot(x,f(x),'--k')
hold on

% Algoritmo di Reinsch
[a1,b1,c1,d1,q] = reinsch(x,y,S,deltay);

s=[];
sval=[];

for i = 1:n-1
    t = x(i):0.00001:x(i+1);
    s=[s, t];
    tval=a1(i)+b1(i)*(t-x(i))+c1(i)*(t-x(i)).^2+d1(i)*(t-x(i)).^3;
    plot(t,tval,'-r');
    hold on
    sval=[sval, tval];
end

lambda = 1/(2*q);

% Algoritmo di Andersen-Chen
p = 2;
gamma = (b-a)^(2*p-1);
[Ut,Vt] = generators1((x-a)./(b-a),p);
[Wt,z] = potrf1(Ut,Vt,n*lambda/gamma);
[ctilde,d] = smoothing_spline_reg(Ut,Wt,z,y);

[phi,zeta]=andersen_chen(a,b,p,n,x);
fhat=@(t) 0;
for k=1:p
    fhat=@(t) fhat(t)+d(k)*phi{k}(t);
end
for j=1:n
    fhat=@(t) fhat(t)+(ctilde(j))/gamma * zeta{j}(t);
end

sval_ac=fhat(s);
plot(s,sval_ac,'--b')
axis([a b -1.5 2])

err=norm(sval-sval_ac)/norm(sval)