function [P,Q,A,G,H,B,d]=quasiseparable(Ut,Wt,c)
% [P,Q,A,G,H,B,d] = quasiseparable(Ut,Wt,c) calcola i generatori
% quasiseparabili di (L*L^T)^{-1}=L^{-T}*L^{-1} dove
% L=tril(Ut'*Wt,-1)+diag(c).

n = size(Ut,2);
[P1,Q1,A1]= inverse_quasisep(Ut,Wt,c);

Q=Q1;
G=Q1';
A=A1;
B=zeros(2*n,2);
for i=1:n
    B(2*i-1:2*i,:)=(A1(2*i-1:2*i,:))';
end

gamma=zeros(2);
P=zeros(n,2);
H=zeros(2,n);
d=zeros(n,1);
for k=n:-1:1
    S=[P1(k,:)', A1(2*k-1:2*k,:)'; 1/c(k), Q1(:,k)'] * [1, zeros(1,2); zeros(2,1), gamma] * [P1(k,:), 1/c(k); A1(2*k-1:2*k,:), Q1(:,k)];
    gamma=S(1:2,1:2);
    H(:,k)=S(1:2,3);
    P(k,:)=S(3,1:2);
    d(k)=S(3,3);
end

end