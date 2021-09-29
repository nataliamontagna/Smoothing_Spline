function [P,Q,A]= inverse_quasisep(Ut,Wt,c)
% [P,Q,A] = inverse_quasisep(Ut,Wt,c) calcola i generatori quasiseparabili
% dell'inversa di L, dove L=tril(Ut'*Wt,-1)+diag(c).

assert(all(size(Ut) == size(Wt)),'Dimension mismatch: Ut and Wt must be of the same size.')
assert(length(c) == size(Ut,2) && isvector(c),'Dimension mismatch: c must be a vector of length size(Ut,2)')

n=size(Ut,2);
P=zeros(n,2);
Q=zeros(2,n);
A=zeros(2*n,2);
for k=1:n
    P(k,:)=-1/c(k) * (Ut(:,k))';
    Q(:,k)=1/c(k) * Wt(:,k);
    A(2*k-1:2*k,:)=eye(2)-1/c(k)* Wt(:,k)*(Ut(:,k))';
end

end