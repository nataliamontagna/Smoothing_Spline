function log_gcv=l_gcv(Ut,Vt,n,gamma,y,v)

[Wt,z]=potrf1(Ut,Vt,n*(10^v)/gamma);
log_gcv=smoothing_spline_gcv(Ut,Wt,z,y);

end