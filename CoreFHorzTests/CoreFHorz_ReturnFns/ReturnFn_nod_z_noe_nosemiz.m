function F=ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension)

F=-Inf;

if agej<Jr
    c=(1+r)*a+w*kappa_j*z-aprime;
else
    c=(1+r)*a+pension-aprime;
end

if c>0
    F=(c^(1-sigma)-1)/(1-sigma);
end


end