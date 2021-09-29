load carbon12alpha

N=length(angle);
n=N-1;
a=0;
b=4.2;

x=angle;
y=counts;

deltay = sqrt(N) * ones(N,1);

S = 150;

plot(x,y,'o')
hold on

[a1,b1,c1,d1,q] = reinsch(x,y,S,deltay);

s=[];
sval=[];

for i = 1:n
    t = x(i):0.001:x(i+1);
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
[Wt,z] = potrf1(Ut,Vt,N*lambda/gamma);
[ctilde,d] = smoothing_spline_reg(Ut,Wt,z,y);

[phi,zeta]=andersen_chen(a,b,p,N,x);
fhat=@(t) 0;
for k=1:p
    fhat=@(t) fhat(t)+d(k)*phi{k}(t);
end
for j=1:N
    fhat=@(t) fhat(t)+(ctilde(j))/gamma * zeta{j}(t);
end

sval_ac=fhat(s);
plot(s,sval_ac,'--b')
axis([a b -50 350])

err=norm(sval-sval_ac)/norm(sval)