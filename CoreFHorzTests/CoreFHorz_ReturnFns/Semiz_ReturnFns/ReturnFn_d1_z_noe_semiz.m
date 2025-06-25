function F=ReturnFn_d1_z_noe_semiz(d1,d2,aprime,a,semiz,z,r,w,kappa_j,sigma,agej,Jr,pension,eta,varphi,uempbenefit,searcheffortcost)

F=-Inf;

if agej<Jr
    c=(1+r)*a+w*kappa_j*z*d1*semiz+uempbenefit*(1-semiz)-aprime;
else
    c=(1+r)*a+pension-aprime;
end

if c>0
    F=(c^(1-sigma)-1)/(1-sigma)+varphi*((1-d1)^(1-eta)-1)/(1-eta)-searcheffortcost*d2;
end


end