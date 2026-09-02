function T=GPTemptationFn_d_noz_e_nosemiz(d,aprime,a,e,lambdaGP,shiftGP,r,w,kappa_j,sigma,agej,Jr,pension)
% Temptation utility: consumption is tempting, v = lambdaGP*u_c(c) + shiftGP (no leisure
% term, so the temptation params deliberately differ from the ReturnFn params: no eta/varphi).
% Feasibility (-Inf) is identical to the matching ReturnFn; the budget is the same.

T=-Inf;

if agej<Jr
    c=(1+r)*a+w*kappa_j*e*d-aprime;
else
    c=(1+r)*a+pension-aprime;
end

if c>0 && d<1
    T=lambdaGP*(c^(1-sigma)-1)/(1-sigma)+shiftGP;
end

end
