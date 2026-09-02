function T=GPTemptationFn_nod1_z_e_semiz(d2,aprime,a,semiz,z,e,lambdaGP,shiftGP,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit)
% Temptation utility: consumption is tempting, v = lambdaGP*u_c(c) + shiftGP (no leisure or
% search-effort terms, so the temptation params deliberately differ from the ReturnFn params).
% Feasibility (-Inf) is identical to the matching ReturnFn; the budget is the same.

T=-Inf;

if agej<Jr
    c=(1+r)*a+w*kappa_j*z*e*semiz+uempbenefit*(1-semiz)-aprime;
else
    c=(1+r)*a+pension-aprime;
end

if c>0
    T=lambdaGP*(c^(1-sigma)-1)/(1-sigma)+shiftGP;
end

end
