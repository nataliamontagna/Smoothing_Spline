function [phi,zeta]=andersen_chen(a,b,p,n,x)

gamma=(b-a)^(2*p-1);
phi=cell(p,1);
for k=1:p
    phi{k}=@(t) (((t-a)./(b-a)).^(k-1))./factorial(k-1);
end

kernel=@(s,t) 0;
for k=0:p-1
    kernel=@(s,t) kernel(s,t) + ((-1)^k)/(factorial(p-1-k)*factorial(p+k)) * (s*t).^(p-1-k) .* (min(s,t)).^(2*k+1);
end

zeta=cell(n,1);
for j=1:n
    zeta{j}=@(t) gamma* kernel((x(j)-a)/(b-a), (t-a)./(b-a) );
end

end