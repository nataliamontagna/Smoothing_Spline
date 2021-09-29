function log_gml=l_gml(Ut,Vt,n,gamma,y,v)

[Wt,z]=potrf1(Ut,Vt,n*(10^v)/gamma);
log_gml=smoothing_spline_par(Ut,Wt,z,y);

end