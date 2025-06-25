function F=ReturnFn_nod1_noz_noe_semiz(d2,aprime,a,semiz,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost)

F=-Inf;

if agej<Jr
    c=(1+r)*a+w*kappa_j*semiz+uempbenefit*(1-semiz)-aprime;
else
    c=(1+r)*a+pension-aprime;
end

if c>0
    F=(c^(1-sigma)-1)/(1-sigma)-searcheffortcost*d2;
end


end