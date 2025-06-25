function F=ReturnFn_d_z_e_nosemiz(d,aprime,a,z,e,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension)

F=-Inf;

if agej<Jr
    c=(1+r)*a+w*kappa_j*z*e*d-aprime;
else
    c=(1+r)*a+pension-aprime;
end

if c>0 && d<1
    F=(c^(1-sigma)-1)/(1-sigma)+varphi*((1-d)^(1-eta)-1)/(1-eta);
end


end