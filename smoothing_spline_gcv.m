function log_gcv=smoothing_spline_gcv(Ut,Wt,z,y)

n=size(Ut,2);
[ctilde,d,Q1,R]=smoothing_spline_reg(Ut,Wt,z,y);

[P2,Q2,A2,G2,H2,B2,f]=quasiseparable(Ut,Wt,z);
traccia = sum(f);

A=zeros(n,2);
A(:,1)=trsv1(Ut,Wt,z,Q1(:,1),'T');
A(:,2)=trsv1(Ut,Wt,z,Q1(:,2),'T');
log_gcv = log10(n) + 2* log10(norm(ctilde,2)) -2* log10(traccia - norm(A,'fro')^2); % A = L'\Q1

end