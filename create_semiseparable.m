function A = create_semiseparable(Ut,Wt,z)
    A = tril(Ut'*Wt,-1)+diag(z);
end