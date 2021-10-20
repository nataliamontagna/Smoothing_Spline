function log_gml=smoothing_spline_gml(Ut,Wt,z,y)

n=size(Ut,2);
[ctilde,d,Q1,R]=smoothing_spline_reg(Ut,Wt,z,y);

log_gml = log10(y'*ctilde) + 2/(n-2) * ( (sum(log10(z))) + sum(log10(abs(diag(R)))) );

end