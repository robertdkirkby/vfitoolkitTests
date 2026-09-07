function output=GPFHorz_CrossTests_d1_semiz_with2A(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)
% The Gul-Pesendorfer with2A semiz cross tests, with d1 (n_d/d_grid passed in are n_d_semiz/d_grid_semiz;
% n_a/a_grid carry the two assets, [101,4]-style with stacked a_grids).
% Test 1: lambdaGP=0 (no temptation) equals the standard semiz solver. Run for all four shock
%         combos (z_noe, noz_e, z_e, noz_noe) at all four tiers (plain, DC, GI, DC+GI), so
%         every GP semiz-2A raw gets a limiting-case anchor.
% Test 2: a semi-exo that is really just a markov, at the BASELINE lambdaGP (temptation fully
%         active on both sides): a SemiExoStateFn that ignores d2 and reproduces pi_z, with
%         semiz_grid=z_grid and uempbenefit=searcheffortcost=0 so the return and temptation
%         fns coincide, must reproduce the (GPU-validated) nosemiz-2A GP solve. At all four
%         tiers. This pins the semiz-2A temptation machinery (per-d2 most-tempting collection,
%         fine-grid refinement under GI) against the known-good nosemiz-2A family.
% Test 3: V_Jplus1 short solve vs full solve (lambdaGP=0.1), at all four tiers (runtime
%         coverage of the GP semiz-2A raws' V_Jplus1 branches).
% Test 4: sign check: temptation weakly lowers welfare, V_GP <= V_standard everywhere.

% Semiz baseline (the employed/not-employed semi-exo state with binary search-effort d2)
n_semiz=vfoptionsbaseline.n_semiz;
semiz_grid=vfoptionsbaseline.semiz_grid;
SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;

n_e=vfoptionsbaseline.n_e;
e_grid=vfoptionsbaseline.e_grid;
pi_e=vfoptionsbaseline.pi_e;

half=N_j/2; % for cross test 3 (N_j=20)

ReturnFn_z=@(d1,d2,a1prime,a2prime,a1,a2,semiz,z,r,w,kappa_j,sigma,agej,Jr,pension,eta,varphi,uempbenefit,searcheffortcost,phi1,phi2) ReturnFn_d1_z_noe_semiz_with2A(d1,d2,a1prime,a2prime,a1,a2,semiz,z,r,w,kappa_j,sigma,agej,Jr,pension,eta,varphi,uempbenefit,searcheffortcost,phi1,phi2);
ReturnFn_e=@(d1,d2,a1prime,a2prime,a1,a2,semiz,e,r,w,kappa_j,sigma,agej,Jr,pension,eta,varphi,uempbenefit,searcheffortcost,phi1,phi2) ReturnFn_d1_noz_e_semiz_with2A(d1,d2,a1prime,a2prime,a1,a2,semiz,e,r,w,kappa_j,sigma,agej,Jr,pension,eta,varphi,uempbenefit,searcheffortcost,phi1,phi2);
ReturnFn_ze=@(d1,d2,a1prime,a2prime,a1,a2,semiz,z,e,r,w,kappa_j,sigma,agej,Jr,pension,eta,varphi,uempbenefit,searcheffortcost,phi1,phi2) ReturnFn_d1_z_e_semiz_with2A(d1,d2,a1prime,a2prime,a1,a2,semiz,z,e,r,w,kappa_j,sigma,agej,Jr,pension,eta,varphi,uempbenefit,searcheffortcost,phi1,phi2);
ReturnFn_none=@(d1,d2,a1prime,a2prime,a1,a2,semiz,r,w,kappa_j,sigma,agej,Jr,pension,eta,varphi,uempbenefit,searcheffortcost,phi1,phi2) ReturnFn_d1_noz_noe_semiz_with2A(d1,d2,a1prime,a2prime,a1,a2,semiz,r,w,kappa_j,sigma,agej,Jr,pension,eta,varphi,uempbenefit,searcheffortcost,phi1,phi2);

TemptationFn_z=@(d1,d2,a1prime,a2prime,a1,a2,semiz,z,lambdaGP,shiftGP,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit) GPTemptationFn_d1_z_noe_semiz_with2A(d1,d2,a1prime,a2prime,a1,a2,semiz,z,lambdaGP,shiftGP,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit);
TemptationFn_e=@(d1,d2,a1prime,a2prime,a1,a2,semiz,e,lambdaGP,shiftGP,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit) GPTemptationFn_d1_noz_e_semiz_with2A(d1,d2,a1prime,a2prime,a1,a2,semiz,e,lambdaGP,shiftGP,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit);
TemptationFn_ze=@(d1,d2,a1prime,a2prime,a1,a2,semiz,z,e,lambdaGP,shiftGP,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit) GPTemptationFn_d1_z_e_semiz_with2A(d1,d2,a1prime,a2prime,a1,a2,semiz,z,e,lambdaGP,shiftGP,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit);
TemptationFn_none=@(d1,d2,a1prime,a2prime,a1,a2,semiz,lambdaGP,shiftGP,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit) GPTemptationFn_d1_noz_noe_semiz_with2A(d1,d2,a1prime,a2prime,a1,a2,semiz,lambdaGP,shiftGP,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit);

% Params with the temptation switched off (used by test 1)
Params0=Params;
Params0.lambdaGP=0;

% The four tiers as vfoptions add-ons, each carrying the semiz fields
vfsemiz=struct();
vfsemiz.n_semiz=n_semiz;
vfsemiz.semiz_grid=semiz_grid;
vfsemiz.SemiExoStateFn=SemiExoStateFn;
vfstd1=vfsemiz;
vfstd2=vfsemiz; vfstd2.divideandconquer=1;
vfstd3=vfsemiz; vfstd3.gridinterplayer=1; vfstd3.ngridinterp=5;
vfstd4=vfstd2; vfstd4.gridinterplayer=1; vfstd4.ngridinterp=5;

%% ==================== Cross test 1: lambdaGP=0 equals the standard semiz solver ====================
% --- combo z_noe ---
vfgp1=vfstd1; vfgp1.exoticpreferences='GulPesendorfer'; vfgp1.temptationFn=TemptationFn_z;
vfgp2=vfstd2; vfgp2.exoticpreferences='GulPesendorfer'; vfgp2.temptationFn=TemptationFn_z;
vfgp3=vfstd3; vfgp3.exoticpreferences='GulPesendorfer'; vfgp3.temptationFn=TemptationFn_z;
vfgp4=vfstd4; vfgp4.exoticpreferences='GulPesendorfer'; vfgp4.temptationFn=TemptationFn_z;
[Vstd1,Pstd1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfstd1);
[Vstd2,Pstd2]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfstd2);
[Vstd3,Pstd3]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfstd3);
[Vstd4,Pstd4]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfstd4);
[V1,P1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params0,DiscountFactorParamNames,[],vfgp1);
fprintf('Cross test 1 (semiz z_noe with2A, plain): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(V1(:)-Vstd1(:))))
fprintf('Cross test 1 (semiz z_noe with2A, plain): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(P1(:)-Pstd1(:))))
[V2,P2]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params0,DiscountFactorParamNames,[],vfgp2);
fprintf('Cross test 1 (semiz z_noe with2A, DC): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(V2(:)-Vstd2(:))))
fprintf('Cross test 1 (semiz z_noe with2A, DC): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(P2(:)-Pstd2(:))))
[V3,P3]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params0,DiscountFactorParamNames,[],vfgp3);
fprintf('Cross test 1 (semiz z_noe with2A, GI): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(V3(:)-Vstd3(:))))
fprintf('Cross test 1 (semiz z_noe with2A, GI): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(P3(:)-Pstd3(:))))
[V4,P4]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params0,DiscountFactorParamNames,[],vfgp4);
fprintf('Cross test 1 (semiz z_noe with2A, DC+GI): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(V4(:)-Vstd4(:))))
fprintf('Cross test 1 (semiz z_noe with2A, DC+GI): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(P4(:)-Pstd4(:))))
clear V2 V3 V4 P2 P3 P4
% (Vstd1..4/Pstd1..4 and V1/P1 for the z_noe combo are kept: tests 3 and 4 reuse them)

% --- combo noz_e ---
vfstd1e=vfsemiz; vfstd1e.n_e=n_e; vfstd1e.e_grid=e_grid; vfstd1e.pi_e=pi_e;
vfstd2e=vfstd1e; vfstd2e.divideandconquer=1;
vfstd3e=vfstd1e; vfstd3e.gridinterplayer=1; vfstd3e.ngridinterp=5;
vfstd4e=vfstd2e; vfstd4e.gridinterplayer=1; vfstd4e.ngridinterp=5;
vfgp1e=vfstd1e; vfgp1e.exoticpreferences='GulPesendorfer'; vfgp1e.temptationFn=TemptationFn_e;
vfgp2e=vfstd2e; vfgp2e.exoticpreferences='GulPesendorfer'; vfgp2e.temptationFn=TemptationFn_e;
vfgp3e=vfstd3e; vfgp3e.exoticpreferences='GulPesendorfer'; vfgp3e.temptationFn=TemptationFn_e;
vfgp4e=vfstd4e; vfgp4e.exoticpreferences='GulPesendorfer'; vfgp4e.temptationFn=TemptationFn_e;
[Vstd,Pstd]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params,DiscountFactorParamNames,[],vfstd1e);
[V1e,P1e]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params0,DiscountFactorParamNames,[],vfgp1e);
fprintf('Cross test 1 (semiz noz_e with2A, plain): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(V1e(:)-Vstd(:))))
fprintf('Cross test 1 (semiz noz_e with2A, plain): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(P1e(:)-Pstd(:))))
[Vstd,Pstd]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params,DiscountFactorParamNames,[],vfstd2e);
[V1e,P1e]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params0,DiscountFactorParamNames,[],vfgp2e);
fprintf('Cross test 1 (semiz noz_e with2A, DC): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(V1e(:)-Vstd(:))))
fprintf('Cross test 1 (semiz noz_e with2A, DC): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(P1e(:)-Pstd(:))))
[Vstd,Pstd]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params,DiscountFactorParamNames,[],vfstd3e);
[V1e,P1e]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params0,DiscountFactorParamNames,[],vfgp3e);
fprintf('Cross test 1 (semiz noz_e with2A, GI): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(V1e(:)-Vstd(:))))
fprintf('Cross test 1 (semiz noz_e with2A, GI): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(P1e(:)-Pstd(:))))
[Vstd,Pstd]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params,DiscountFactorParamNames,[],vfstd4e);
[V1e,P1e]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params0,DiscountFactorParamNames,[],vfgp4e);
fprintf('Cross test 1 (semiz noz_e with2A, DC+GI): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(V1e(:)-Vstd(:))))
fprintf('Cross test 1 (semiz noz_e with2A, DC+GI): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(P1e(:)-Pstd(:))))
clear V1e P1e Vstd Pstd

% --- combo z_e ---
vfgp1ze=vfstd1e; vfgp1ze.exoticpreferences='GulPesendorfer'; vfgp1ze.temptationFn=TemptationFn_ze;
vfgp2ze=vfstd2e; vfgp2ze.exoticpreferences='GulPesendorfer'; vfgp2ze.temptationFn=TemptationFn_ze;
vfgp3ze=vfstd3e; vfgp3ze.exoticpreferences='GulPesendorfer'; vfgp3ze.temptationFn=TemptationFn_ze;
vfgp4ze=vfstd4e; vfgp4ze.exoticpreferences='GulPesendorfer'; vfgp4ze.temptationFn=TemptationFn_ze;
[Vstd,Pstd]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_ze,Params,DiscountFactorParamNames,[],vfstd1e);
[V1ze,P1ze]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_ze,Params0,DiscountFactorParamNames,[],vfgp1ze);
fprintf('Cross test 1 (semiz z_e with2A, plain): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(V1ze(:)-Vstd(:))))
fprintf('Cross test 1 (semiz z_e with2A, plain): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(P1ze(:)-Pstd(:))))
[Vstd,Pstd]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_ze,Params,DiscountFactorParamNames,[],vfstd2e);
[V1ze,P1ze]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_ze,Params0,DiscountFactorParamNames,[],vfgp2ze);
fprintf('Cross test 1 (semiz z_e with2A, DC): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(V1ze(:)-Vstd(:))))
fprintf('Cross test 1 (semiz z_e with2A, DC): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(P1ze(:)-Pstd(:))))
[Vstd,Pstd]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_ze,Params,DiscountFactorParamNames,[],vfstd3e);
[V1ze,P1ze]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_ze,Params0,DiscountFactorParamNames,[],vfgp3ze);
fprintf('Cross test 1 (semiz z_e with2A, GI): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(V1ze(:)-Vstd(:))))
fprintf('Cross test 1 (semiz z_e with2A, GI): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(P1ze(:)-Pstd(:))))
[Vstd,Pstd]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_ze,Params,DiscountFactorParamNames,[],vfstd4e);
[V1ze,P1ze]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_ze,Params0,DiscountFactorParamNames,[],vfgp4ze);
fprintf('Cross test 1 (semiz z_e with2A, DC+GI): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(V1ze(:)-Vstd(:))))
fprintf('Cross test 1 (semiz z_e with2A, DC+GI): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(P1ze(:)-Pstd(:))))
clear V1ze P1ze Vstd Pstd

% --- combo noz_noe ---
vfgp1n=vfstd1; vfgp1n.exoticpreferences='GulPesendorfer'; vfgp1n.temptationFn=TemptationFn_none;
vfgp2n=vfstd2; vfgp2n.exoticpreferences='GulPesendorfer'; vfgp2n.temptationFn=TemptationFn_none;
vfgp3n=vfstd3; vfgp3n.exoticpreferences='GulPesendorfer'; vfgp3n.temptationFn=TemptationFn_none;
vfgp4n=vfstd4; vfgp4n.exoticpreferences='GulPesendorfer'; vfgp4n.temptationFn=TemptationFn_none;
[Vstd,Pstd]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_none,Params,DiscountFactorParamNames,[],vfstd1);
[V1n,P1n]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_none,Params0,DiscountFactorParamNames,[],vfgp1n);
fprintf('Cross test 1 (semiz noz_noe with2A, plain): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(V1n(:)-Vstd(:))))
fprintf('Cross test 1 (semiz noz_noe with2A, plain): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(P1n(:)-Pstd(:))))
[Vstd,Pstd]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_none,Params,DiscountFactorParamNames,[],vfstd2);
[V1n,P1n]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_none,Params0,DiscountFactorParamNames,[],vfgp2n);
fprintf('Cross test 1 (semiz noz_noe with2A, DC): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(V1n(:)-Vstd(:))))
fprintf('Cross test 1 (semiz noz_noe with2A, DC): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(P1n(:)-Pstd(:))))
[Vstd,Pstd]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_none,Params,DiscountFactorParamNames,[],vfstd3);
[V1n,P1n]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_none,Params0,DiscountFactorParamNames,[],vfgp3n);
fprintf('Cross test 1 (semiz noz_noe with2A, GI): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(V1n(:)-Vstd(:))))
fprintf('Cross test 1 (semiz noz_noe with2A, GI): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(P1n(:)-Pstd(:))))
[Vstd,Pstd]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_none,Params,DiscountFactorParamNames,[],vfstd4);
[V1n,P1n]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_none,Params0,DiscountFactorParamNames,[],vfgp4n);
fprintf('Cross test 1 (semiz noz_noe with2A, DC+GI): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(V1n(:)-Vstd(:))))
fprintf('Cross test 1 (semiz noz_noe with2A, DC+GI): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(P1n(:)-Pstd(:))))
clear V1n P1n Vstd Pstd

%% ==================== Cross test 2: a semi-exo that is really just a markov (baseline lambdaGP) ====================
% The SemiExoStateFn ignores d2 and reproduces pi_z exactly, with semiz_grid=z_grid and
% uempbenefit=searcheffortcost=0 so the return and temptation fns coincide with the nosemiz
% z-model. Temptation is ACTIVE (baseline lambdaGP), so this pins the GP semiz machinery
% against the GPU-validated nosemiz GP family. The semiz Policy carries an extra (irrelevant)
% d2 row, which is dropped before comparing.
n_dz=n_d(1); % just the d1 part (the z-side model has no d2)
d_gridz=d_grid(1:n_d(1));
Params2=Params;
Params2.uempbenefit=0;
Params2.searcheffortcost=0;
n_zJAM=2;
z_gridJAM=[0.6;1.4];
pi_zJAM=[1-Params.probfindjob, Params.probfindjob; Params.problosejob, 1-Params.problosejob];
Params2.z1=z_gridJAM(1);
Params2.z2=z_gridJAM(2);
SemiExoStateFn_JustAMarkov=@(n,nprime,dsemiz,probfindjob,problosejob,z1,z2) CoreFHorzSetup_SemiExoStateFn_JustAMarkov(n,nprime,dsemiz,probfindjob,problosejob,z1,z2);
% The nosemiz z-side (uses the GPU-validated nosemiz-2A GP family; with d, z plays the semiz role)
ReturnFn_zside=@(d,a1prime,a2prime,a1,a2,z,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension,phi1,phi2) ReturnFn_d_z_noe_nosemiz_with2A(d,a1prime,a2prime,a1,a2,z,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension,phi1,phi2);
TemptationFn_zside=@(d,a1prime,a2prime,a1,a2,z,lambdaGP,shiftGP,r,w,kappa_j,sigma,agej,Jr,pension) GPTemptationFn_d_z_noe_nosemiz_with2A(d,a1prime,a2prime,a1,a2,z,lambdaGP,shiftGP,r,w,kappa_j,sigma,agej,Jr,pension);
% The semiz side: no z, the semi-exo state IS the markov (uses ReturnFn_none/TemptationFn_none with uempbenefit=0, searcheffortcost=0)
vfJAM=struct();
vfJAM.n_semiz=n_zJAM;
vfJAM.semiz_grid=z_gridJAM;
vfJAM.SemiExoStateFn=SemiExoStateFn_JustAMarkov;
for tier_c=1:4
    vfz=struct();
    vfz.exoticpreferences='GulPesendorfer';
    vfz.temptationFn=TemptationFn_zside;
    vfs=vfJAM;
    vfs.exoticpreferences='GulPesendorfer';
    vfs.temptationFn=TemptationFn_none;
    if tier_c==1
        tiername='plain';
    elseif tier_c==2
        tiername='DC';
        vfz.divideandconquer=1; vfs.divideandconquer=1;
    elseif tier_c==3
        tiername='GI';
        vfz.gridinterplayer=1; vfz.ngridinterp=5; vfs.gridinterplayer=1; vfs.ngridinterp=5;
    elseif tier_c==4
        tiername='DC+GI';
        vfz.divideandconquer=1; vfz.gridinterplayer=1; vfz.ngridinterp=5;
        vfs.divideandconquer=1; vfs.gridinterplayer=1; vfs.ngridinterp=5;
    end
    [Vz,Pz]=ValueFnIter_Case1_FHorz(n_dz,n_a,n_zJAM,N_j,d_gridz,a_grid,z_gridJAM,pi_zJAM,ReturnFn_zside,Params2,DiscountFactorParamNames,[],vfz);
    [Vs,Ps]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_none,Params2,DiscountFactorParamNames,[],vfs);
    fprintf('Cross test 2 (semiz just-a-markov with2A, %s): V vs nosemiz GP, this should be zero: %.3e \n',tiername,max(abs(Vs(:)-Vz(:))))
    % Drop the d2 row of the semiz Policy (row 2 here: rows are (d1,d2,a1prime,a2prime,...)), then compare
    Psshort=Ps([1,3:size(Ps,1)],:,:,:);
    fprintf('Cross test 2 (semiz just-a-markov with2A, %s): Policy (without d2 row) vs nosemiz GP, this should be zero: %.3e \n',tiername,max(abs(Psshort(:)-Pz(:))))
end
clear Vz Vs Pz Ps Psshort

%% ==================== Cross test 3: V_Jplus1 short solve vs full solve (baseline lambdaGP) ====================
% (the z_noe combo; the full-horizon GP solves at each tier vs a short solve fed V_Jplus1)
Njs=half;
Paramsjs=Params;
Paramsjs.agej=Params.agej(1:Njs);
Paramsjs.kappa_j=Params.kappa_j(1:Njs);
% plain tier (V1/P1 are the lambda=0 runs; re-solve at baseline lambda for this test)
[V5a,P5a]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfgp1);
vfgp1c=vfgp1;
vfgp1c.V_Jplus1=V5a(:,:,:,:,half+1);
[V5c,P5c]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Paramsjs,DiscountFactorParamNames,[],vfgp1c);
temp=reshape(V5a,[],N_j); temp=temp(:,1:half)-reshape(V5c,[],Njs);
fprintf('Cross test 3 (semiz z_noe with2A, plain): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
temp=reshape(P5a,[],N_j); temp=temp(:,1:half)-reshape(P5c,[],Njs);
fprintf('Cross test 3 (semiz z_noe with2A, plain): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
% DC tier
[V5a,P5a]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfgp2);
vfgp2c=vfgp2;
vfgp2c.V_Jplus1=V5a(:,:,:,:,half+1);
[V5c,P5c]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Paramsjs,DiscountFactorParamNames,[],vfgp2c);
temp=reshape(V5a,[],N_j); temp=temp(:,1:half)-reshape(V5c,[],Njs);
fprintf('Cross test 3 (semiz z_noe with2A, DC): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
temp=reshape(P5a,[],N_j); temp=temp(:,1:half)-reshape(P5c,[],Njs);
fprintf('Cross test 3 (semiz z_noe with2A, DC): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
% GI tier
[V5a,P5a]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfgp3);
vfgp3c=vfgp3;
vfgp3c.V_Jplus1=V5a(:,:,:,:,half+1);
[V5c,P5c]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Paramsjs,DiscountFactorParamNames,[],vfgp3c);
temp=reshape(V5a,[],N_j); temp=temp(:,1:half)-reshape(V5c,[],Njs);
fprintf('Cross test 3 (semiz z_noe with2A, GI): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
temp=reshape(P5a,[],N_j); temp=temp(:,1:half)-reshape(P5c,[],Njs);
fprintf('Cross test 3 (semiz z_noe with2A, GI): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
% DC+GI tier
[V5a,P5a]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfgp4);
vfgp4c=vfgp4;
vfgp4c.V_Jplus1=V5a(:,:,:,:,half+1);
[V5c,P5c]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Paramsjs,DiscountFactorParamNames,[],vfgp4c);
temp=reshape(V5a,[],N_j); temp=temp(:,1:half)-reshape(V5c,[],Njs);
fprintf('Cross test 3 (semiz z_noe with2A, DC+GI): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
temp=reshape(P5a,[],N_j); temp=temp(:,1:half)-reshape(P5c,[],Njs);
fprintf('Cross test 3 (semiz z_noe with2A, DC+GI): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
clear V5a V5c P5a P5c temp

%% ==================== Cross test 4: temptation weakly lowers welfare ====================
% Self-control cost is nonnegative, so V_GP <= V_standard everywhere (baseline lambdaGP)
[Vgp,~]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfgp1);
temp=max(Vgp(:)-Vstd1(:));
fprintf('Cross test 4 (semiz z_noe with2A, plain): max(V_GP - V_standard), this should be weakly negative: %.3e \n',temp)

%%
output=struct(); % Not currently used for anything. Maybe will do so later.

end
