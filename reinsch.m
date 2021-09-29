function [a,b,c,d,p,cont]=reinsch(x,y,S,deltay)
% [a,b,c,d,p,cont]=reinsch(x,y,S,deltay) calcola i coefficienti della 
% spline che risolve il problema della smoothing spline regression, il
% parametro lagrangiano p e il numero di iterazioni dell'algoritmo cont

n=length(x)-1;
h=zeros(n,1);
h=x(2:n+1)-x(1:n);

D=spdiags(deltay,0,n+1,n+1);

T_diag=2*(h(1:n-1)+h(2:n))/3;
T_sdiag=h(2:n-1)/3;
B=[[T_sdiag; 0], T_diag, [0; T_sdiag]];
T=spdiags(B,-1:1,n-1,n-1);

Q_diag=1./h(1:n-1);
Q_sdiag=-1./h(1:n-1)-1./h(2:n);
Q_ssdiag=1./h(2:n);
Q=spdiags([Q_ssdiag, Q_sdiag, Q_diag], -2:0, n+1, n-1);

p=0;
M=(Q')*(D^2)*Q;
R=chol(M);

u=M\(Q'*y);
v=D*Q*u;
e=v'*v;
cont=0;
while e>S
    f=u'*T*u;
    w=R'\(T*u);
    g=w'*w;
    p=p+(e-sqrt(S*e))/(f-p*g);
    M=(Q')*(D^2)*Q+p*T;
    R=chol(M);
    u=M\(Q'*y);
    v=D*Q*u;
    e=v'*v;
    cont=cont+1;
end

a=y-D*v;
c=zeros(n+1,1);
c(2:n)=p*u;

d=(c(2:n+1)-c(1:n))./(3*h(1:n));

b=(a(2:n+1)-a(1:n))./h(1:n) - c(1:n).*h(1:n) - d(1:n).*(h(1:n)).^2;

a=a(1:n);
c=c(1:n);

end