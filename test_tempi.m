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

% Algoritmo di Andersen-Chen
p = 2;
gamma = (b-a)^(2*p-1);
lambda = 1.e-9;
[Ut,Vt] = generators1((x-a)./(b-a),p);
[Wt,z] = potrf1(Ut,Vt,n*lambda/gamma);
y = f(x) + sigma * randn(n,1);

tic;
[ctilde,d] = smoothing_spline_reg(Ut,Wt,z,y);
yhat=y-n*lambda/gamma * ctilde;
t1=toc

% Algoritmo di Reinsch
m = n-1;

h=zeros(m,1);
h=x(2:m+1)-x(1:m);

D=spdiags(deltay,0,n,n);

T_diag=2*(h(1:m-1)+h(2:m))/3;
T_sdiag=h(2:m-1)/3;
B=[[T_sdiag; 0], T_diag, [0; T_sdiag]];
T=spdiags(B,-1:1,m-1,m-1);

Q_diag=1./h(1:m-1);
Q_sdiag=-1./h(1:m-1)-1./h(2:m);
Q_ssdiag=1./h(2:m);
Q=spdiags([Q_ssdiag, Q_sdiag, Q_diag], -2:0, m+1, m-1);

M = Q'*D^2*Q + 1/(2*lambda)*T;
w = M\(Q'*y);
S = norm(D*Q*w)^2;

tic;
[a1,b1,c1,d1,q,cont] = reinsch(x,y,S,deltay);
t2=toc
t3=t2/cont