function output=AmbRiskyAsset_nod1_noz_noe_nosemiz_withA1(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c)
% AmbiguityAversion mirror of CoreFHorzRiskyAsset_nod1_noz_noe_nosemiz_withA1. V/Policy/ValueFnFromPolicy only
% (exotic preferences change nothing downstream of Policy). Ambiguity is over pi_u (mandatory: the risky
% return distribution is ambiguous, not known risk).

% Setup vfoptions
vfoptions=struct();
% Riskyasset
vfoptions.riskyasset=1;
vfoptions.refine_d=[0,1,1]; % no d1, d2 (riskyshare), d3 (savings)
vfoptions.aprimeFn=vfoptionsbaseline.aprimeFn;
vfoptions.n_u=vfoptionsbaseline.n_u;
vfoptions.u_grid=vfoptionsbaseline.u_grid;
vfoptions.pi_u=vfoptionsbaseline.pi_u;
% Do the current setup
n_z=0;
z_grid=[];
pi_z=[];

ReturnFn=@(savings,a1prime,a1,a2,r,w,kappa_j,sigma,r_a1,agej,Jr,pension) ReturnFn_nod1_noz_noe_nosemiz_withA1(savings,a1prime,a1,a2,r,w,kappa_j,sigma,r_a1,agej,Jr,pension);

%% Ambiguity Aversion (u is treated as ambiguity, not risk: ambiguity_pi_u is mandatory)
vfoptions.exoticpreferences='AmbiguityAversion';
vfoptions.n_ambiguity=vfoptionsbaseline.n_ambiguity;
vfoptions.ambiguity_pi_u=vfoptionsbaseline.ambiguity_pi_u;
% The regular pi_u/pi_z/pi_e (used for the agent distribution) equal prior 1, so no consistency warning

%%
vfoptions1=vfoptions;
[V1,Policy1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);

% Solve with divide-and-conquer, should give same answer
vfoptions2=vfoptions;
vfoptions2.divideandconquer=1;
[V2,Policy2]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);
fprintf('Divide-and-conquer, this should be zero: %.3e \n',max(abs(V1(:)-V2(:))))
fprintf('Divide-and-conquer, this should be zero: %.3e \n',max(abs(Policy1(:)-Policy2(:))))

% V from Policy
V1fromPolicy=ValueFnFromPolicy_FHorz(Policy1,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions1);
fprintf('ValueFnFromPolicy, this should be zero: %.3e \n',max(abs(V1fromPolicy(:)-V1(:))))

% lowmemory
vfoptions1.lowmemory=1;
[V1B,Policy1B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
fprintf('lowmemory=1, this should be zero: %.3e \n',max(abs(V1(:)-V1B(:))))
fprintf('lowmemory=1, this should be zero: %.3e \n',max(abs(Policy1(:)-Policy1B(:))))
vfoptions1.lowmemory=0;

vfoptions2.lowmemory=1;
[V2B,Policy2B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);
fprintf('lowmemory=1 (with DC), this should be zero: %.3e \n',max(abs(V2(:)-V2B(:))))
fprintf('lowmemory=1 (with DC), this should be zero: %.3e \n',max(abs(Policy2(:)-Policy2B(:))))
vfoptions2.lowmemory=0;

clear V1 V2 V1B V2B Policy1 Policy2 Policy1B Policy2B V1fromPolicy

%% Solve with grid-interpolation
vfoptions3=vfoptions;
vfoptions3.gridinterplayer=1;
vfoptions3.ngridinterp=5;
[V3,Policy3]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);

% Solve with divide-and-conquer, should give same answer
vfoptions4=vfoptions;
vfoptions4.divideandconquer=1;
vfoptions4.gridinterplayer=1;
vfoptions4.ngridinterp=5;
[V4,Policy4]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions4);

fprintf('Divide-and-conquer (with Grid Interp Layer), this should be zero: %.3e \n',max(abs(V3(:)-V4(:))))
fprintf('Divide-and-conquer (with Grid Interp Layer), this should be zero: %.3e \n',max(abs(Policy3(:)-Policy4(:))))

% V from Policy
V3fromPolicy=ValueFnFromPolicy_FHorz(Policy3,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions3);
fprintf('ValueFnFromPolicy (GI), this should be zero: %.3e \n',max(abs(V3fromPolicy(:)-V3(:))))

% lowmemory
vfoptions3.lowmemory=1;
[V3B,Policy3B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);
fprintf('lowmemory=1 (with GI), this should be zero: %.3e \n',max(abs(V3(:)-V3B(:))))
fprintf('lowmemory=1 (with GI), this should be zero: %.3e \n',max(abs(Policy3(:)-Policy3B(:))))
vfoptions3.lowmemory=0;

vfoptions4.lowmemory=1;
[V4B,Policy4B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions4);
fprintf('lowmemory=1 (with DC+GI), this should be zero: %.3e \n',max(abs(V4(:)-V4B(:))))
fprintf('lowmemory=1 (with DC+GI), this should be zero: %.3e \n',max(abs(Policy4(:)-Policy4B(:))))
vfoptions4.lowmemory=0;

%%
output=struct(); % Not currently used for anything. Maybe will do so later.

end
