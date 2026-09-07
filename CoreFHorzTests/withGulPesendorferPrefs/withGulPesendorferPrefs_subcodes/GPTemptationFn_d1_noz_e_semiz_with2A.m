function T=GPTemptationFn_d1_noz_e_semiz_with2A(d1,d2,a1prime,a2prime,a1,a2,semiz,e,lambdaGP,shiftGP,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit)
% Temptation utility: consumption is tempting, v = lambdaGP*u_c(c) + shiftGP (no leisure,
% search-effort or a2-preference terms, so the temptation params deliberately differ from
% the ReturnFn params). Feasibility (-Inf) is identical to the matching ReturnFn; the budget
% is the same.

T=-Inf;

if agej<Jr
    c=(1+r)*a1+w*kappa_j*d1*semiz*e+uempbenefit*(1-semiz)-a1prime+a2-a2prime;
else
    c=(1+r)*a1+pension-a1prime+a2-a2prime;
end

if c>0 && d1<1
    T=lambdaGP*(c^(1-sigma)-1)/(1-sigma)+shiftGP;
end

end
