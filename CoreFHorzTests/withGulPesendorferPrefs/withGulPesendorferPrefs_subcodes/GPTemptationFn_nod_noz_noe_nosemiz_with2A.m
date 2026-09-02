function T=GPTemptationFn_nod_noz_noe_nosemiz_with2A(a1prime,a2prime,a1,a2,lambdaGP,shiftGP,r,w,kappa_j,sigma,agej,Jr,pension)
% Temptation utility: consumption is tempting, v = lambdaGP*u_c(c) + shiftGP (no leisure or
% a2-preference terms, so the temptation params deliberately differ from the ReturnFn params:
% no eta/varphi/phi1/phi2). Feasibility (-Inf) is identical to the matching ReturnFn; the
% budget is the same.

T=-Inf;

if agej<Jr
    c=(1+r)*a1+w*kappa_j-a1prime+a2-a2prime;
else
    c=(1+r)*a1+pension-a1prime+a2-a2prime;
end

if c>0
    T=lambdaGP*(c^(1-sigma)-1)/(1-sigma)+shiftGP;
end

end
