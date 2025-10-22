function F=ReturnFn_TwoEndo_d1_noz_noe_nosemiz(d1,a1prime,a2prime,a1,a2,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension)

F=-Inf;

if agej<Jr
    c=(1+r)*a1+w*kappa_j*d1*a2prime*a2-a1prime;
else
    c=(1+r)*a1+pension-a1prime;
end

if c>0 && d1<1
    F=(c^(1-sigma)-1)/(1-sigma)+varphi*((1-d1)^(1-eta)-1)/(1-eta);
end


end