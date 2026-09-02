function T=GPTemptationFn_nod_noz_noe_nosemiz(aprime,a,lambdaGP,shiftGP,r,w,kappa_j,sigma,agej,Jr,pension)
% Temptation utility: consumption is tempting, v = lambdaGP*u_c(c) + shiftGP.
% Feasibility (-Inf) is identical to the matching ReturnFn; the budget is the same.

T=-Inf;

if agej<Jr
    c=(1+r)*a+w*kappa_j-aprime;
else
    c=(1+r)*a+pension-aprime;
end

if c>0
    T=lambdaGP*(c^(1-sigma)-1)/(1-sigma)+shiftGP;
end

end
