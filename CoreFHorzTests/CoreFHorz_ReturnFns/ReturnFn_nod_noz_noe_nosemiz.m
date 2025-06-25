function F=ReturnFn_nod_noz_noe_nosemiz(aprime,a,r,w,kappa_j,sigma,agej,Jr,pension)

F=-Inf;

if agej<Jr
    c=(1+r)*a+w*kappa_j-aprime;
else
    c=(1+r)*a+pension-aprime;
end

if c>0
    F=(c^(1-sigma)-1)/(1-sigma);
end


end