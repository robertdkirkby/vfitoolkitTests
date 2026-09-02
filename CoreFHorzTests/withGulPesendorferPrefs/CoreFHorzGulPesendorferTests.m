% Implement lots of tests of the core VFI Toolkit FHorz commands, under GUL-PESENDORFER
% preferences (temptation and self-control: a temptation fn v alongside the return fn u, with
%   V_j = max_{d,a'} [ u + v + beta*E V_{j+1} ] - max_{d,a'} v
% so Policy maximizes the tempted objective and V nets off the self-control cost)
% with/without d
% with/without z and/or e (Gul-Pesendorfer is fine with no shocks at all, so all 8 combos)
% with/without divide-and-conquer
% with/without grid interpolation
% with/without low memory (where appropriate; the z&e models also get lowmemory=2)
%
% This is the Gul-Pesendorfer mirror of CoreFHorzTests.m. The temptation fn is declared as
% vfoptions.temptationFn, with the same input signature as the ReturnFn (its parameters may
% differ; here they deliberately do: the temptation is over consumption only). Baseline
% temptation: v = lambdaGP*u_c(c) + shiftGP with lambdaGP=0.1, shiftGP=0 (a la Krusell,
% Kuruscu & Smith); shiftGP exists purely so the cross tests can exercise the
% constant-cancellation and shift-invariance properties.
%
% TEST-FIRST STATE (2026-09-02, see GulPesendorfer_FHorz_proposal.md in the toolkit repo):
% the whole bank is test-first, written against the planned full 1A GulPesendorfer family.
% Against the pre-existing toolkit: the plain tier has 7 of 8 raws (the noz-with-e raw's
% dispatch line is commented out, so the d_noz_e case errors); the DC/GI/DC+GI flags are
% silently IGNORED by the old GP dispatcher (so old-toolkit DC checks would pass trivially and
% GI solves return a wrong-shaped Policy); and ValueFnFromPolicy has no GulPesendorfer
% branch (its generic path omits the self-control cost, so those checks would be nonzero).
% All of these are fixed by the accompanying toolkit wave; run this bank against it.
% No with2A section yet (deliberately deferred, unlike the QH/Ambiguity banks).

%% Diary of the command window output (figures are saved into the same folder as they are created)
if ~exist('../TestOutput','dir')
    mkdir('../TestOutput')
end
if exist('../TestOutput/CoreFHorzGulPesendorferTestsdiary.txt','file')
    delete('../TestOutput/CoreFHorzGulPesendorferTestsdiary.txt') % otherwise diary just appends to the previous run
end
diary ../TestOutput/CoreFHorzGulPesendorferTestsdiary.txt

addpath('../CoreFHorzTests_Setup/')
addpath('../CoreFHorz_ReturnFns/')

addpath('./withGulPesendorferPrefs_subcodes/')
addpath('./withGulPesendorferPrefs_subcodes/CrossTests/')

% Setup so that use the same d,a,z,e,semiz in all the models that use them
CoreFHorz_setup

%% Gul-Pesendorfer: the temptation parameters
% The temptation fns live in withGulPesendorferPrefs_subcodes/ (GPTemptationFn_*), one per
% case signature; each subcode declares its own vfoptions.temptationFn (mirroring how the
% ReturnFn is declared per-subcode)
Params.lambdaGP=0.1; % temptation strength: v = lambdaGP*u_c(c) + shiftGP
Params.shiftGP=0; % constant shift of the temptation utility; should never change anything

% vfoptions.exoticpreferences='GulPesendorfer' is set inside each subcode

%% The only functions worth testing are the value fn ones, as after you have Policy everything else is anyway ignoring the temptation

%% Gul-Pesendorfer without declaring vfoptions.temptationFn is deliberately an error
vfoptionstemp.exoticpreferences='GulPesendorfer';
ReturnFn_none=@(aprime,a,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_noz_noe_nosemiz(aprime,a,r,w,kappa_j,sigma,agej,Jr,pension);
try
    [Vtemp,Policytemp]=ValueFnIter_Case1_FHorz(0,n_a,0,N_j,[],a_grid,[],[],ReturnFn_none,Params,DiscountFactorParamNames,[],vfoptionstemp);
    fprintf('GulPesendorfer without temptationFn: FAIL, this should error and did not \n')
catch ME
    if contains(ME.message,'temptationFn')
        fprintf('GulPesendorfer without temptationFn errors as intended :) \n')
    else
        fprintf('GulPesendorfer without temptationFn: FAIL, errored but with the wrong message: %s \n',ME.message)
    end
end
clear vfoptionstemp

%% Gul-Pesendorfer with semi-exogenous states is deliberately an error (semiz tier not implemented)
vfoptionstemp.exoticpreferences='GulPesendorfer';
vfoptionstemp.temptationFn=@(aprime,a,lambdaGP,shiftGP,r,w,kappa_j,sigma,agej,Jr,pension) GPTemptationFn_nod_noz_noe_nosemiz(aprime,a,lambdaGP,shiftGP,r,w,kappa_j,sigma,agej,Jr,pension);
vfoptionstemp.n_semiz=vfoptionsbaseline.n_semiz;
vfoptionstemp.semiz_grid=vfoptionsbaseline.semiz_grid;
vfoptionstemp.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;
try
    [Vtemp,Policytemp]=ValueFnIter_Case1_FHorz(n_d2_semiz,n_a,0,N_j,d2_grid_semiz,a_grid,[],[],ReturnFn_none,Params,DiscountFactorParamNames,[],vfoptionstemp);
    fprintf('GulPesendorfer with semiz: FAIL, this should error and did not \n')
catch ME
    if contains(ME.message,'semi')
        fprintf('GulPesendorfer with semiz errors as intended :) \n')
    else
        fprintf('GulPesendorfer with semiz: FAIL, errored but with the wrong message: %s \n',ME.message)
    end
end
clear vfoptionstemp

%% without d, with z, without e
figure_c=1;
output=GPFHorz_nod_z_noe_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzGulPesendorferTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with d, with z, without e
figure_c=2;
output=GPFHorz_d_z_noe_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzGulPesendorferTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% without d, without z, with e
figure_c=3;
output=GPFHorz_nod_noz_e_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzGulPesendorferTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with d, without z, with e
figure_c=4;
output=GPFHorz_d_noz_e_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzGulPesendorferTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% without d, with z, with e
figure_c=5;
output=GPFHorz_nod_z_e_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzGulPesendorferTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with d, with z, with e
figure_c=6;
output=GPFHorz_d_z_e_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzGulPesendorferTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% without d, without z, without e (Gul-Pesendorfer works with no shocks at all)
figure_c=7;
output=GPFHorz_nod_noz_noe_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzGulPesendorferTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with d, without z, without e
figure_c=8;
output=GPFHorz_d_noz_noe_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzGulPesendorferTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% The cross tests (see the comments at the top of the cross-test subcodes for what they cover)
output=GPFHorz_CrossTests_nod_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

output=GPFHorz_CrossTests_d_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

%% Done
diary off
