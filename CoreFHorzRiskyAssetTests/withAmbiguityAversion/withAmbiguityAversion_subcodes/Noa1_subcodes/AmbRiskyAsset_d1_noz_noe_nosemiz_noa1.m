function output=AmbRiskyAsset_d1_noz_noe_nosemiz_noa1(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c)
% AmbiguityAversion mirror of CoreFHorzRiskyAsset_d1_noz_noe_nosemiz_noa1. V/Policy/ValueFnFromPolicy only
% (exotic preferences change nothing downstream of Policy). Ambiguity is over pi_u (mandatory: the risky
% return distribution is ambiguous, not known risk).
% noa1 riskyasset has no divide-and-conquer and no grid interpolation (no a1 to refine), so only the
% base method appears here; the four-method legs are in the WithA1 subcodes.

% Setup vfoptions
vfoptions=struct();
% Riskyasset
vfoptions.riskyasset=1;
vfoptions.refine_d=[1,1,1]; % d1 (h), d2 (riskyshare), d3 (savings)
vfoptions.aprimeFn=vfoptionsbaseline.aprimeFn;
vfoptions.n_u=vfoptionsbaseline.n_u;
vfoptions.u_grid=vfoptionsbaseline.u_grid;
vfoptions.pi_u=vfoptionsbaseline.pi_u;
% Do the current setup
n_z=0;
z_grid=[];
pi_z=[];

ReturnFn=@(h,savings,a,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension) ReturnFn_d1_noz_noe_nosemiz(h,savings,a,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension);

%% Ambiguity Aversion (u is treated as ambiguity, not risk: ambiguity_pi_u is mandatory)
vfoptions.exoticpreferences='AmbiguityAversion';
vfoptions.n_ambiguity=vfoptionsbaseline.n_ambiguity;
vfoptions.ambiguity_pi_u=vfoptionsbaseline.ambiguity_pi_u;
% The regular pi_u/pi_z/pi_e (used for the agent distribution) equal prior 1, so no consistency warning

%%
vfoptions1=vfoptions;
[V1,Policy1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);

% V from Policy
V1fromPolicy=ValueFnFromPolicy_FHorz(Policy1,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions1);
fprintf('ValueFnFromPolicy, this should be zero: %.3e \n',max(abs(V1fromPolicy(:)-V1(:))))

% lowmemory
vfoptions1.lowmemory=1;
[V1B,Policy1B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
fprintf('lowmemory=1, this should be zero: %.3e \n',max(abs(V1(:)-V1B(:))))
fprintf('lowmemory=1, this should be zero: %.3e \n',max(abs(Policy1(:)-Policy1B(:))))
vfoptions1.lowmemory=0;

%%
output=struct(); % Not currently used for anything. Maybe will do so later.

end
