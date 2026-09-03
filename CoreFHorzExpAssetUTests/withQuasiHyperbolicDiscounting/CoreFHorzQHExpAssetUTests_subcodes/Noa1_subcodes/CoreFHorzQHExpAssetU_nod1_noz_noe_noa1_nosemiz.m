function output=CoreFHorzQHExpAssetU_nod1_noz_noe_noa1_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c)
% Quasi-hyperbolic experienceassetu, noa1: no d1, no z, no e.
% The experience asset a2 is the only endogenous state. n_a is scalar (n_a_justexpasset);
% a_grid is the a2_grid. n_a_big/a_grid_big are unused. n_d=n_d2; d_grid=d2_grid
%
% Methods: BASE only -- with no a1 there is nothing for divide-and-conquer or the grid
% interpolation layer to operate on (same reason the baseline noa1 subcodes have no
% DC/GI/DC+GI blocks). Runs Naive then Sophisticated, each at lowmemory=0 only (there
% are no shocks, so no higher lowmemory level exists) and with a ValueFnFromPolicy
% oracle, then the exponential-discounting cross-tests.
% shocks: {} (none) -> valid lowmemory {0}.
%
% TEST-FIRST: the toolkit currently has NO quasi-hyperbolic support for experienceassetu at
% all (no QH+ExpAssetU raws, no dispatcher branch), so this errors at the first ValueFnIter
% call. That is expected: this test is written ahead of the toolkit code.

% Setup vfoptions and simoptions
vfoptions=struct();
simoptions=struct();
n_z=0;
z_grid=[];
pi_z=[];

ReturnFn=@(d2,a,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod1_noz_noe_noa1_nosemiz(d2,a,r,w,kappa_j,sigma,agej,Jr,pension);

% Experience asset u
vfoptions.experienceassetu=1;
simoptions.experienceassetu=1;
vfoptions.aprimeFn=vfoptionsbaseline.aprimeFn;
simoptions.aprimeFn=vfoptions.aprimeFn;
simoptions.d_grid=d_grid;
simoptions.a_grid=a_grid;
vfoptions.n_u=vfoptionsbaseline.n_u;
vfoptions.u_grid=vfoptionsbaseline.u_grid;
vfoptions.pi_u=vfoptionsbaseline.pi_u;
simoptions.n_u=vfoptions.n_u;
simoptions.u_grid=vfoptions.u_grid;
simoptions.pi_u=vfoptions.pi_u;

%% Quasi-Hyperbolic Discounting
vfoptions.exoticpreferences='QuasiHyperbolic';
vfoptions.QHadditionaldiscount=vfoptionsbaseline.QHadditionaldiscount;

%% Naive
vfoptions.quasi_hyperbolic='Naive';
vfoptions1=vfoptions;
[V1,Policy1,V1alt,Policy1alt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);

% 0 shocks (no z, no e, no semiz) -> only lowmemory=0 is valid; no lowmemory>0 check.

vfoptions1.Policyalt=Policy1alt; % Naive QH: ValueFnFromPolicy reconstructs V from the exponential-discounter argmax (the 4th output of the Naive solve)
[V1fromPolicy,V1altfromPolicy]=ValueFnFromPolicy_FHorz(Policy1,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions1);
fprintf('Naive ValueFnFromPolicy, this should be zero: %.3e \n',max(abs(V1fromPolicy(:)-V1(:))))
fprintf('Naive ValueFnFromPolicy (Valt), this should be zero: %.3e \n',max(abs(V1altfromPolicy(:)-V1alt(:))))

%% Sophisticated
vfoptions.quasi_hyperbolic='Sophisticated';
vfoptions1=vfoptions;
[V1s,Policy1s,V1salt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);

% 0 shocks (no z, no e, no semiz) -> only lowmemory=0 is valid; no lowmemory>0 check.

vfoptions1.Policyalt=[]; % Sophisticated does not need Policyalt
[V1sfromPolicy,V1saltfromPolicy]=ValueFnFromPolicy_FHorz(Policy1s,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions1);
fprintf('Sophisticated ValueFnFromPolicy, this should be zero: %.3e \n',max(abs(V1sfromPolicy(:)-V1s(:))))
fprintf('Sophisticated ValueFnFromPolicy (Valt), this should be zero: %.3e \n',max(abs(V1saltfromPolicy(:)-V1salt(:))))


%% V_Jplus1: use Valt of period jstar as the terminal value function of a shorter model
% As in the other test banks: solve the model, then solve a shorter model that runs only periods
% 1,...,jstar-1 with Njs=jstar-1 and the age-dependent parameters trimmed to length Njs, and check
% we get the same answer for periods 1,...,jstar-1.
% IMPORTANT (quasi-hyperbolic specific): the continuation value that enters the Bellman equation is
% the STANDARD-discounted one -- V_std when Naive, Vunderbar when Sophisticated -- which is the
% third output Valt. It is NOT V (which is Vtilde/Vhat, the quasi-hyperbolic-discounted object).
% So it is Valt that gets fed back in as vfoptions.V_Jplus1.
% Note: this block sits before the "Versus exponential discounting" section on purpose -- that
% section sets Params.beta0=1, which would collapse QH to exponential and stop these tests from
% exercising any of the quasi-hyperbolic machinery.
% Note: vfoptions1..4 were last redefined in the Sophisticated section above, so each block below
% sets quasi_hyperbolic explicitly rather than relying on what they currently hold.
% Note: u is iid and vfoptions.pi_u carries no age dimension, so it needs no trimming; u is
% integrated out inside the a2prime transition and so adds no dimension to V or Policy.
% Note: mewj is age-dependent, but is only used for the agent distribution, which is not computed
% here, so it is left alone.
% This tier has no a1 for divide-and-conquer or the grid interp layer to operate on, so it
% defines only vfoptions1. Run at jstar=round(3*N_j/4) and again at jstar=N_j, so the
% retirement periods and the terminal V_Jplus1 branch are covered too.

%% V_Jplus1, Naive
for jstar=[round(3*N_j/4),N_j]
    Njs=jstar-1; % the shorter model runs periods 1,...,jstar-1
    Paramsjs=Params;
    Paramsjs.agej=Params.agej(1:Njs);
    Paramsjs.kappa_j=Params.kappa_j(1:Njs);
    vfoptionsjs=vfoptions1;
    vfoptionsjs.exoticpreferences='QuasiHyperbolic';
    vfoptionsjs.quasi_hyperbolic='Naive';
    [Vbase,Policybase,Valtbase,Policyaltbase]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsjs);
    vfoptionsjs.V_Jplus1=Valtbase(:,jstar); % Valt, not V
    Vbase=Vbase(:,1:Njs);
    Policybase=Policybase(:,:,1:Njs);
    Valtbase=Valtbase(:,1:Njs);
    Policyaltbase=Policyaltbase(:,:,1:Njs);
    [Vshort,Policyshort,Valtshort,Policyaltshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
    fprintf('V_Jplus1 (jstar=%i, Naive), this should be zero: %.3e \n',jstar,max(abs(Vbase(:)-Vshort(:))))
    fprintf('V_Jplus1 (jstar=%i, Naive), this should be zero: %.3e \n',jstar,max(abs(Policybase(:)-Policyshort(:))))
    fprintf('V_Jplus1 (jstar=%i, Naive, Valt), this should be zero: %.3e \n',jstar,max(abs(Valtbase(:)-Valtshort(:))))
    fprintf('V_Jplus1 (jstar=%i, Naive, Policyalt), this should be zero: %.3e \n',jstar,max(abs(Policyaltbase(:)-Policyaltshort(:))))
end

%% V_Jplus1, Sophisticated
for jstar=[round(3*N_j/4),N_j]
    Njs=jstar-1; % the shorter model runs periods 1,...,jstar-1
    Paramsjs=Params;
    Paramsjs.agej=Params.agej(1:Njs);
    Paramsjs.kappa_j=Params.kappa_j(1:Njs);
    vfoptionsjs=vfoptions1;
    vfoptionsjs.exoticpreferences='QuasiHyperbolic';
    vfoptionsjs.quasi_hyperbolic='Sophisticated';
    [Vbase,Policybase,Valtbase]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsjs);
    vfoptionsjs.V_Jplus1=Valtbase(:,jstar); % Valt, not V
    Vbase=Vbase(:,1:Njs);
    Policybase=Policybase(:,:,1:Njs);
    Valtbase=Valtbase(:,1:Njs);
    [Vshort,Policyshort,Valtshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
    fprintf('V_Jplus1 (jstar=%i, Sophisticated), this should be zero: %.3e \n',jstar,max(abs(Vbase(:)-Vshort(:))))
    fprintf('V_Jplus1 (jstar=%i, Sophisticated), this should be zero: %.3e \n',jstar,max(abs(Policybase(:)-Policyshort(:))))
    fprintf('V_Jplus1 (jstar=%i, Sophisticated, Valt), this should be zero: %.3e \n',jstar,max(abs(Valtbase(:)-Valtshort(:))))
end

clear Vbase Policybase Valtbase Policyaltbase Vshort Policyshort Valtshort Policyaltshort

%% Versus exponential discounting
% (i) at the actual QH beta0, Naive's continuation value equals the exponential value function
vfoptionsE=vfoptions; vfoptionsE.exoticpreferences='None';
[Vexp,Policyexp]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsE);
vfoptionsN=vfoptions; vfoptionsN.quasi_hyperbolic='Naive';
[VnA,PolicynA,VnAalt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsN);
fprintf('(i) Naive continuation value == exponential (beta0=%g), should be zero: %.3e \n',Params.beta0,max(abs(VnAalt(:)-Vexp(:))))

% (ii) with beta0=1, Naive main value AND continuation value both equal exponential
% (iii) with beta0=1, Sophisticated main value AND continuation value both equal exponential
beta0_store=Params.beta0; Params.beta0=1;
[Vexp1,~]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsE);
[VnA1,~,VnA1alt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsN);
vfoptionsS=vfoptions; vfoptionsS.quasi_hyperbolic='Sophisticated';
[VsS1,~,VsS1alt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsS);
fprintf('(ii) Naive V == exponential (beta0=1), should be zero: %.3e \n',max(abs(VnA1(:)-Vexp1(:))))
fprintf('(ii) Naive Valt == exponential (beta0=1), should be zero: %.3e \n',max(abs(VnA1alt(:)-Vexp1(:))))
fprintf('(iii) Sophisticated V == exponential (beta0=1), should be zero: %.3e \n',max(abs(VsS1(:)-Vexp1(:))))
fprintf('(iii) Sophisticated Valt == exponential (beta0=1), should be zero: %.3e \n',max(abs(VsS1alt(:)-Vexp1(:))))
Params.beta0=beta0_store;

%%
output=struct(); % Not currently used for anything. Maybe will do so later.

end
