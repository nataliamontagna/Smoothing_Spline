function [ctilde,d]=stability_and(Ut,Wt,z,y)
n=size(Ut,2);
L=create_semiseparable(Ut,Wt,z);
k=cond(L)

B=zeros(n,2);
B(:,1)=trsv1(Ut,Wt,z,Ut(2,:)','N');
B(:,2)=trsv1(Ut,Wt,z,Ut(1,:)','N');

[v,w,Q1,R]=qr_householder(B);
R1=R(1:2,:);
ctilde=trsv1(Ut,Wt,z,y,'N');
ctilde1=L\y;
ctilde2=vpa(L,32)\vpa(y,32);

err_indietro=double(norm(vpa(y,32)-vpa(L,32)*vpa(ctilde,32))/norm(vpa(y,32)))
err_indietro1=double(norm(vpa(y,32)-vpa(L,32)*vpa(ctilde1,32))/norm(vpa(y,32)))
err_avanti=double(norm(vpa(ctilde,32)-ctilde2)/norm(ctilde2))
err_avanti1=double(norm(vpa(ctilde1,32)-ctilde2)/norm(ctilde2))
t=ctilde-2*(v'*ctilde)*v; ctilde=t-2*(w'*t)*w;

d=R1\(ctilde(1:2));
ctilde(1:2)=0;
t1=ctilde-2*(w'*ctilde)*w; s=t1-2*(v'*t1)*v;

k1=cond(L')
ctilde=trsv1(Ut,Wt,z,s,'T');
ctilde1=L'\s;
ctilde2=vpa(L',32)\vpa(s);
err_indietro_=double(norm(vpa(s,32)-vpa(L',32)*vpa(ctilde,32))/norm(vpa(s,32)))
err_indietro_1=double(norm(vpa(s,32)-vpa(L',32)*vpa(ctilde1,32))/norm(vpa(s,32)))
err_avanti_=double(norm(vpa(ctilde,32)-ctilde2)/norm(ctilde2))
err_avanti_1=double(norm(vpa(ctilde1,32)-ctilde2)/norm(ctilde2))

end