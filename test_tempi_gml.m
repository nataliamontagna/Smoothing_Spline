n = 500;
a = 0;
b = 1;

sigma = 0.1;
x = zeros(n,1);
for i=1:n
    x(i)=(i-1)/(n-1);
end
f = @(t) cos(2*pi*t) + 0.3*sin(10*pi*t);
deltay = sqrt(n) * ones(n,1);

p=2;
gamma=(b-a)^(2*p-1);
[Ut,Vt]=generators1((x-a)./(b-a),p);

t=zeros(10,1);
for k=1:10

y = f(x) + sigma * randn(n,1);

tic;
q=fminbnd(@(v) l_gml(Ut,Vt,n,gamma,y,v), -10, 10);
lambda=10^q;
t(k)=toc;

end
tmean=mean(t)