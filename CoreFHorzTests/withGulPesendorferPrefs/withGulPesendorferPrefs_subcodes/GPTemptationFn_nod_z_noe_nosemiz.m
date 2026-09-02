function T=GPTemptationFn_nod_z_noe_nosemiz(aprime,a,z,lambdaGP,shiftGP,r,w,kappa_j,sigma,agej,Jr,pension)
% Temptation utility: consumption is tempting, v = lambdaGP*u_c(c) + shiftGP.
% Feasibility (-Inf) is identical to the matching ReturnFn; the budget is the same.
% shiftGP is a constant shift of the temptation utility, which should never change
% anything (it cancels against the most-tempting term); the cross tests exploit this.

T=-Inf;

if agej<Jr
    c=(1+r)*a+w*kappa_j*z-aprime;
else
    c=(1+r)*a+pension-aprime;
end

if c>0
    T=lambdaGP*(c^(1-sigma)-1)/(1-sigma)+shiftGP;
end

end
