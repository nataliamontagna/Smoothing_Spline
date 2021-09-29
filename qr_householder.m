function [v,w,Q1,R]=qr_householder(B)

n=size(B,1);
x=B(:,1);
alpha=norm(x,2);
e1=zeros(n,1);
e1(1)=1;
u=x-alpha*e1;
v=u./norm(u,2);

b=B(:,1)-2*(v'*B(:,1)).*v;
b11=b(1);
b1=b(2:end,:);

c=B(:,2)-2*(v'*B(:,2)).*v;
b12=c(1);
b2=c(2:end,:);

alpha1=norm(b2,2);
e1=e1(1:n-1);
u1=b2-alpha1*e1;
v1=u1./norm(u1,2);
w=[0; v1];

I=eye(n);
e1=I(:,1);
e2=I(:,2);
Q1=zeros(n,2);
Q1(:,1)=e1+(4*(v'*w)*w(1)-2*v(1))*v-2*w(1)*w;
Q1(:,2)=e2+(4*(v'*w)*w(2)-2*v(2))*v-2*w(2)*w;

R=zeros(n,2);
R(1,:)=[b11 b12];
R(2:end,1)=b1-2*(v1'*b1)*v1;
R(2:end,2)=b2-2*(v1'*b2)*v1;
end