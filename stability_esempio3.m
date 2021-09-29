load carbon12alpha

n=length(angle);
a=0;
b=4.2;

x=angle;
y=counts;

p=2;
gamma=(b-a)^(2*p-1);
[Ut,Vt]=generators1((x-a)./(b-a),p);

lambda=1.e-1;

[Wt,z]=potrf1(Ut,Vt,n*lambda/gamma);

[ctilde,d]=stability_and(Ut,Wt,z,y);