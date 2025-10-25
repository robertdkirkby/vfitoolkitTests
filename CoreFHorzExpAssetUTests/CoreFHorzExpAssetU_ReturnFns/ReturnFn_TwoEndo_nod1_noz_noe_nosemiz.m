function F=ReturnFn_TwoEndo_nod1_noz_noe_nosemiz(a1prime,a2prime,a1,a2,r,w,kappa_j,sigma,agej,Jr,pension)

F=-Inf;

if agej<Jr
    c=(1+r)*a1+w*kappa_j*a2prime*a2-a1prime;
else
    c=(1+r)*a1+pension-a1prime;
end

if c>0
    F=(c^(1-sigma)-1)/(1-sigma);
end


end