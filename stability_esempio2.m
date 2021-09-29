Y=load('mcyc.dat');
tobs=Y(:,1); yobs=Y(:,2);
[t,y,w]=splav(tobs,yobs);

n=length(t);
a=0;
b=60;

p=2;
gamma=(b-a)^(2*p-1);
[Ut,Vt]=generators1((t-a)./(b-a),p);

lambda=1.e-10;

[Wt,z]=potrf1(Ut,Vt,n*lambda/gamma);

[ctilde,d]=stability_and(Ut,Wt,z,y);