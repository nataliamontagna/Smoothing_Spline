function [ctilde,d,Q1,R]=smoothing_spline_reg(Ut,Wt,z,y)
% [ctilde,d,Q1,R]=smoothing_spline_reg(Ut,Wt,z,y) risolve il sistema alla
% base del metodo strutturato di Andersen e Chen

n=size(Ut,2);
B=zeros(n,2);
B(:,1)=trsv1(Ut,Wt,z,Ut(2,:)','N');
B(:,2)=trsv1(Ut,Wt,z,Ut(1,:)','N');

[v,w,Q1,R]=qr_householder(B);
R1=R(1:2,:);
ctilde=trsv1(Ut,Wt,z,y,'N');
t=ctilde-2*(v'*ctilde)*v; ctilde=t-2*(w'*t)*w; %ctilde=Q'*ctilde;
d=R1\(ctilde(1:2));
ctilde(1:2)=[0;0];
t1=ctilde-2*(w'*ctilde)*w; s=t1-2*(v'*t1)*v; %s=Q*ctilde;
ctilde=trsv1(Ut,Wt,z,s,'T');

end