function T=GPTemptationFn_nod1_noz_noe_semiz_with2A(d2,a1prime,a2prime,a1,a2,semiz,lambdaGP,shiftGP,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit)
% Temptation utility: consumption is tempting, v = lambdaGP*u_c(c) + shiftGP (no leisure,
% search-effort or a2-preference terms, so the temptation params deliberately differ from
% the ReturnFn params). Feasibility (-Inf) is identical to the matching ReturnFn; the budget
% is the same.

T=-Inf;

if agej<Jr
    c=(1+r)*a1+w*kappa_j*semiz+uempbenefit*(1-semiz)-a1prime+a2-a2prime;
else
    c=(1+r)*a1+pension-a1prime+a2-a2prime;
end

if c>0
    T=lambdaGP*(c^(1-sigma)-1)/(1-sigma)+shiftGP;
end

end
