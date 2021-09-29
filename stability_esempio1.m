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

lambda=1.e-13;

[Wt,z]=potrf1(Ut,Vt,n*lambda/gamma);

[ctilde,d]=stability_and(Ut,Wt,z,y);