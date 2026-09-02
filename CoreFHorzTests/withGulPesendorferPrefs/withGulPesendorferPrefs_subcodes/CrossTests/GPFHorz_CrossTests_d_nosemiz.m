function output=GPFHorz_CrossTests_d_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)
% The Gul-Pesendorfer cross tests, with a decision variable d. Same test list as the nod
% cross tests (see GPFHorz_CrossTests_nod_nosemiz.m for the full comments) except that the
% hand-rolled brute force maximizes jointly over (d,aprime) — the most-tempting term is a max
% over the full JOINT choice set — and there is no hand-rolled GI here (test 4b is nod-only).
% Note: the temptation fns deliberately have different parameters to the ReturnFns (no
% eta/varphi), so these tests also cover TemptationFnParamNames differing from ReturnFnParamNames.
% Test 1: lambdaGP=0 equals standard; all four shock combos at all four tiers.
% Test 2: constant temptation (lambdaGP=0, shiftGP=0.7) equals standard.
% Test 3: shift invariance: v and v+0.7 give identical V and Policy.
% Test 4: hand-rolled brute force (joint (d,aprime) max, full-choice-set most-tempting term)
%         vs the toolkit at the plain and DC tiers.
% Test 5: age-dependent lambdaGP via V_Jplus1, at all four tiers.
% Test 6: sign check: V_GP <= V_standard everywhere.

ReturnFn_z=@(d,aprime,a,z,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension) ReturnFn_d_z_noe_nosemiz(d,aprime,a,z,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension);
ReturnFn_e=@(d,aprime,a,e,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension) ReturnFn_d_noz_e_nosemiz(d,aprime,a,e,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension);
ReturnFn_ze=@(d,aprime,a,z,e,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension) ReturnFn_d_z_e_nosemiz(d,aprime,a,z,e,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension);
ReturnFn_none=@(d,aprime,a,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension) ReturnFn_d_noz_noe_nosemiz(d,aprime,a,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension);

TemptationFn_z=@(d,aprime,a,z,lambdaGP,shiftGP,r,w,kappa_j,sigma,agej,Jr,pension) GPTemptationFn_d_z_noe_nosemiz(d,aprime,a,z,lambdaGP,shiftGP,r,w,kappa_j,sigma,agej,Jr,pension);
TemptationFn_e=@(d,aprime,a,e,lambdaGP,shiftGP,r,w,kappa_j,sigma,agej,Jr,pension) GPTemptationFn_d_noz_e_nosemiz(d,aprime,a,e,lambdaGP,shiftGP,r,w,kappa_j,sigma,agej,Jr,pension);
TemptationFn_ze=@(d,aprime,a,z,e,lambdaGP,shiftGP,r,w,kappa_j,sigma,agej,Jr,pension) GPTemptationFn_d_z_e_nosemiz(d,aprime,a,z,e,lambdaGP,shiftGP,r,w,kappa_j,sigma,agej,Jr,pension);
TemptationFn_none=@(d,aprime,a,lambdaGP,shiftGP,r,w,kappa_j,sigma,agej,Jr,pension) GPTemptationFn_d_noz_noe_nosemiz(d,aprime,a,lambdaGP,shiftGP,r,w,kappa_j,sigma,agej,Jr,pension);

n_e=vfoptionsbaseline.n_e;
e_grid=vfoptionsbaseline.e_grid;
pi_e=vfoptionsbaseline.pi_e;

half=N_j/2; % for cross test 5 (N_j=20)

% Params with the temptation switched off (used by tests 1 and 2)
Params0=Params;
Params0.lambdaGP=0;

%% ==================== Cross test 1: lambdaGP=0 equals standard ====================
% --- combo z_noe ---
vfstd1=struct();
vfstd2=struct(); vfstd2.divideandconquer=1;
vfstd3=struct(); vfstd3.gridinterplayer=1; vfstd3.ngridinterp=5;
vfstd4=struct(); vfstd4.divideandconquer=1; vfstd4.gridinterplayer=1; vfstd4.ngridinterp=5;
vfgp1=vfstd1; vfgp1.exoticpreferences='GulPesendorfer'; vfgp1.temptationFn=TemptationFn_z;
vfgp2=vfstd2; vfgp2.exoticpreferences='GulPesendorfer'; vfgp2.temptationFn=TemptationFn_z;
vfgp3=vfstd3; vfgp3.exoticpreferences='GulPesendorfer'; vfgp3.temptationFn=TemptationFn_z;
vfgp4=vfstd4; vfgp4.exoticpreferences='GulPesendorfer'; vfgp4.temptationFn=TemptationFn_z;
[Vstd1,Pstd1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfstd1);
[Vstd2,Pstd2]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfstd2);
[Vstd3,Pstd3]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfstd3);
[Vstd4,Pstd4]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfstd4);
[V1,P1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params0,DiscountFactorParamNames,[],vfgp1);
fprintf('Cross test 1 (z_noe, plain): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(V1(:)-Vstd1(:))))
fprintf('Cross test 1 (z_noe, plain): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(P1(:)-Pstd1(:))))
[V2,P2]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params0,DiscountFactorParamNames,[],vfgp2);
fprintf('Cross test 1 (z_noe, DC): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(V2(:)-Vstd2(:))))
fprintf('Cross test 1 (z_noe, DC): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(P2(:)-Pstd2(:))))
[V3,P3]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params0,DiscountFactorParamNames,[],vfgp3);
fprintf('Cross test 1 (z_noe, GI): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(V3(:)-Vstd3(:))))
fprintf('Cross test 1 (z_noe, GI): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(P3(:)-Pstd3(:))))
[V4,P4]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params0,DiscountFactorParamNames,[],vfgp4);
fprintf('Cross test 1 (z_noe, DC+GI): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(V4(:)-Vstd4(:))))
fprintf('Cross test 1 (z_noe, DC+GI): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(P4(:)-Pstd4(:))))
clear V1 V2 V3 V4 P1 P2 P3 P4
% (Vstd1..4/Pstd1..4 for the z_noe combo are kept: tests 2, 5 and 6 reuse them)

% --- combo noz_e ---
vfstd1e=struct(); vfstd1e.n_e=n_e; vfstd1e.e_grid=e_grid; vfstd1e.pi_e=pi_e;
vfstd2e=vfstd1e; vfstd2e.divideandconquer=1;
vfstd3e=vfstd1e; vfstd3e.gridinterplayer=1; vfstd3e.ngridinterp=5;
vfstd4e=vfstd2e; vfstd4e.gridinterplayer=1; vfstd4e.ngridinterp=5;
vfgp1e=vfstd1e; vfgp1e.exoticpreferences='GulPesendorfer'; vfgp1e.temptationFn=TemptationFn_e;
vfgp2e=vfstd2e; vfgp2e.exoticpreferences='GulPesendorfer'; vfgp2e.temptationFn=TemptationFn_e;
vfgp3e=vfstd3e; vfgp3e.exoticpreferences='GulPesendorfer'; vfgp3e.temptationFn=TemptationFn_e;
vfgp4e=vfstd4e; vfgp4e.exoticpreferences='GulPesendorfer'; vfgp4e.temptationFn=TemptationFn_e;
[Vstd,Pstd]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params,DiscountFactorParamNames,[],vfstd1e);
[V1,P1]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params0,DiscountFactorParamNames,[],vfgp1e);
fprintf('Cross test 1 (noz_e, plain): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(V1(:)-Vstd(:))))
fprintf('Cross test 1 (noz_e, plain): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(P1(:)-Pstd(:))))
[Vstd,Pstd]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params,DiscountFactorParamNames,[],vfstd2e);
[V1,P1]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params0,DiscountFactorParamNames,[],vfgp2e);
fprintf('Cross test 1 (noz_e, DC): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(V1(:)-Vstd(:))))
fprintf('Cross test 1 (noz_e, DC): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(P1(:)-Pstd(:))))
[Vstd,Pstd]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params,DiscountFactorParamNames,[],vfstd3e);
[V1,P1]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params0,DiscountFactorParamNames,[],vfgp3e);
fprintf('Cross test 1 (noz_e, GI): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(V1(:)-Vstd(:))))
fprintf('Cross test 1 (noz_e, GI): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(P1(:)-Pstd(:))))
[Vstd,Pstd]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params,DiscountFactorParamNames,[],vfstd4e);
[V1,P1]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params0,DiscountFactorParamNames,[],vfgp4e);
fprintf('Cross test 1 (noz_e, DC+GI): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(V1(:)-Vstd(:))))
fprintf('Cross test 1 (noz_e, DC+GI): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(P1(:)-Pstd(:))))
clear V1 P1 Vstd Pstd

% --- combo z_e ---
vfstd1ze=struct(); vfstd1ze.n_e=n_e; vfstd1ze.e_grid=e_grid; vfstd1ze.pi_e=pi_e;
vfstd2ze=vfstd1ze; vfstd2ze.divideandconquer=1;
vfstd3ze=vfstd1ze; vfstd3ze.gridinterplayer=1; vfstd3ze.ngridinterp=5;
vfstd4ze=vfstd2ze; vfstd4ze.gridinterplayer=1; vfstd4ze.ngridinterp=5;
vfgp1ze=vfstd1ze; vfgp1ze.exoticpreferences='GulPesendorfer'; vfgp1ze.temptationFn=TemptationFn_ze;
vfgp2ze=vfstd2ze; vfgp2ze.exoticpreferences='GulPesendorfer'; vfgp2ze.temptationFn=TemptationFn_ze;
vfgp3ze=vfstd3ze; vfgp3ze.exoticpreferences='GulPesendorfer'; vfgp3ze.temptationFn=TemptationFn_ze;
vfgp4ze=vfstd4ze; vfgp4ze.exoticpreferences='GulPesendorfer'; vfgp4ze.temptationFn=TemptationFn_ze;
[Vstd,Pstd]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_ze,Params,DiscountFactorParamNames,[],vfstd1ze);
[V1,P1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_ze,Params0,DiscountFactorParamNames,[],vfgp1ze);
fprintf('Cross test 1 (z_e, plain): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(V1(:)-Vstd(:))))
fprintf('Cross test 1 (z_e, plain): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(P1(:)-Pstd(:))))
[Vstd,Pstd]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_ze,Params,DiscountFactorParamNames,[],vfstd2ze);
[V1,P1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_ze,Params0,DiscountFactorParamNames,[],vfgp2ze);
fprintf('Cross test 1 (z_e, DC): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(V1(:)-Vstd(:))))
fprintf('Cross test 1 (z_e, DC): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(P1(:)-Pstd(:))))
[Vstd,Pstd]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_ze,Params,DiscountFactorParamNames,[],vfstd3ze);
[V1,P1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_ze,Params0,DiscountFactorParamNames,[],vfgp3ze);
fprintf('Cross test 1 (z_e, GI): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(V1(:)-Vstd(:))))
fprintf('Cross test 1 (z_e, GI): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(P1(:)-Pstd(:))))
[Vstd,Pstd]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_ze,Params,DiscountFactorParamNames,[],vfstd4ze);
[V1,P1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_ze,Params0,DiscountFactorParamNames,[],vfgp4ze);
fprintf('Cross test 1 (z_e, DC+GI): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(V1(:)-Vstd(:))))
fprintf('Cross test 1 (z_e, DC+GI): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(P1(:)-Pstd(:))))
clear V1 P1 Vstd Pstd

% --- combo noz_noe ---
vfstd1n=struct();
vfstd2n=struct(); vfstd2n.divideandconquer=1;
vfstd3n=struct(); vfstd3n.gridinterplayer=1; vfstd3n.ngridinterp=5;
vfstd4n=struct(); vfstd4n.divideandconquer=1; vfstd4n.gridinterplayer=1; vfstd4n.ngridinterp=5;
vfgp1n=vfstd1n; vfgp1n.exoticpreferences='GulPesendorfer'; vfgp1n.temptationFn=TemptationFn_none;
vfgp2n=vfstd2n; vfgp2n.exoticpreferences='GulPesendorfer'; vfgp2n.temptationFn=TemptationFn_none;
vfgp3n=vfstd3n; vfgp3n.exoticpreferences='GulPesendorfer'; vfgp3n.temptationFn=TemptationFn_none;
vfgp4n=vfstd4n; vfgp4n.exoticpreferences='GulPesendorfer'; vfgp4n.temptationFn=TemptationFn_none;
[Vstd,Pstd]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_none,Params,DiscountFactorParamNames,[],vfstd1n);
[V1,P1]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_none,Params0,DiscountFactorParamNames,[],vfgp1n);
fprintf('Cross test 1 (noz_noe, plain): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(V1(:)-Vstd(:))))
fprintf('Cross test 1 (noz_noe, plain): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(P1(:)-Pstd(:))))
[Vstd,Pstd]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_none,Params,DiscountFactorParamNames,[],vfstd2n);
[V1,P1]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_none,Params0,DiscountFactorParamNames,[],vfgp2n);
fprintf('Cross test 1 (noz_noe, DC): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(V1(:)-Vstd(:))))
fprintf('Cross test 1 (noz_noe, DC): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(P1(:)-Pstd(:))))
[Vstd,Pstd]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_none,Params,DiscountFactorParamNames,[],vfstd3n);
[V1,P1]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_none,Params0,DiscountFactorParamNames,[],vfgp3n);
fprintf('Cross test 1 (noz_noe, GI): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(V1(:)-Vstd(:))))
fprintf('Cross test 1 (noz_noe, GI): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(P1(:)-Pstd(:))))
[Vstd,Pstd]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_none,Params,DiscountFactorParamNames,[],vfstd4n);
[V1,P1]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_none,Params0,DiscountFactorParamNames,[],vfgp4n);
fprintf('Cross test 1 (noz_noe, DC+GI): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(V1(:)-Vstd(:))))
fprintf('Cross test 1 (noz_noe, DC+GI): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(P1(:)-Pstd(:))))
clear V1 P1 Vstd Pstd

%% ==================== Cross test 2: constant temptation equals standard ====================
% v = 0.7 on the feasible set (lambdaGP=0, shiftGP=0.7): the constant cancels exactly against
% the most-tempting term, so this equals standard preferences DESPITE the temptation machinery
% being fully engaged (nonzero MostTempting). Catches a dropped or double-counted subtraction.
Params2=Params;
Params2.lambdaGP=0;
Params2.shiftGP=0.7;
[V1,P1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params2,DiscountFactorParamNames,[],vfgp1);
fprintf('Cross test 2 (z_noe, plain): constant temptation equals standard, this should be zero: %.3e \n',max(abs(V1(:)-Vstd1(:))))
fprintf('Cross test 2 (z_noe, plain): constant temptation equals standard, this should be zero: %.3e \n',max(abs(P1(:)-Pstd1(:))))
[V2,P2]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params2,DiscountFactorParamNames,[],vfgp2);
fprintf('Cross test 2 (z_noe, DC): constant temptation equals standard, this should be zero: %.3e \n',max(abs(V2(:)-Vstd2(:))))
fprintf('Cross test 2 (z_noe, DC): constant temptation equals standard, this should be zero: %.3e \n',max(abs(P2(:)-Pstd2(:))))
[V3,P3]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params2,DiscountFactorParamNames,[],vfgp3);
fprintf('Cross test 2 (z_noe, GI): constant temptation equals standard, this should be zero: %.3e \n',max(abs(V3(:)-Vstd3(:))))
fprintf('Cross test 2 (z_noe, GI): constant temptation equals standard, this should be zero: %.3e \n',max(abs(P3(:)-Pstd3(:))))
[V4,P4]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params2,DiscountFactorParamNames,[],vfgp4);
fprintf('Cross test 2 (z_noe, DC+GI): constant temptation equals standard, this should be zero: %.3e \n',max(abs(V4(:)-Vstd4(:))))
fprintf('Cross test 2 (z_noe, DC+GI): constant temptation equals standard, this should be zero: %.3e \n',max(abs(P4(:)-Pstd4(:))))
clear V1 V2 V3 V4 P1 P2 P3 P4

%% ==================== Cross test 3: shift invariance ====================
% v and v+0.7 (lambdaGP=0.1, shiftGP 0 vs 0.7) give identical V and Policy
Params3a=Params; % lambdaGP=0.1, shiftGP=0 (the baseline)
Params3b=Params;
Params3b.shiftGP=0.7;
[V1a,P1a]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params3a,DiscountFactorParamNames,[],vfgp1);
[V1b,P1b]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params3b,DiscountFactorParamNames,[],vfgp1);
fprintf('Cross test 3 (z_noe, plain): shift invariance, this should be zero: %.3e \n',max(abs(V1a(:)-V1b(:))))
fprintf('Cross test 3 (z_noe, plain): shift invariance, this should be zero: %.3e \n',max(abs(P1a(:)-P1b(:))))
[V2a,P2a]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params3a,DiscountFactorParamNames,[],vfgp2);
[V2b,P2b]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params3b,DiscountFactorParamNames,[],vfgp2);
fprintf('Cross test 3 (z_noe, DC): shift invariance, this should be zero: %.3e \n',max(abs(V2a(:)-V2b(:))))
fprintf('Cross test 3 (z_noe, DC): shift invariance, this should be zero: %.3e \n',max(abs(P2a(:)-P2b(:))))
[V3a,P3a]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params3a,DiscountFactorParamNames,[],vfgp3);
[V3b,P3b]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params3b,DiscountFactorParamNames,[],vfgp3);
fprintf('Cross test 3 (z_noe, GI): shift invariance, this should be zero: %.3e \n',max(abs(V3a(:)-V3b(:))))
fprintf('Cross test 3 (z_noe, GI): shift invariance, this should be zero: %.3e \n',max(abs(P3a(:)-P3b(:))))
[V4a,P4a]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params3a,DiscountFactorParamNames,[],vfgp4);
[V4b,P4b]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params3b,DiscountFactorParamNames,[],vfgp4);
fprintf('Cross test 3 (z_noe, DC+GI): shift invariance, this should be zero: %.3e \n',max(abs(V4a(:)-V4b(:))))
fprintf('Cross test 3 (z_noe, DC+GI): shift invariance, this should be zero: %.3e \n',max(abs(P4a(:)-P4b(:))))
clear V2a V2b V3a V3b V4a V4b P2a P2b P3a P3b P4a P4b
% (V1a/P1a, the plain-tier GP solve at the baseline lambdaGP=0.1, is kept for test 6)

%% ==================== Cross test 4: hand-rolled brute force (plain and DC) ====================
% Small grid, plain MATLAB loops on the CPU, calling the very same ReturnFn/TemptationFn anons.
% The max is JOINT over (d,aprime) (rows ordered with d fastest, matching the toolkit's kron
% ordering), and the most-tempting term is the max over the full joint choice set. Operation
% order matches the raws: max((U+T)+beta*EV) with the most-tempting term subtracted AFTER the
% max. Toolkit Policy rows (d index, aprime index) are recombined into the joint index
% d + n_d*(aprime-1) for the comparison.
n_a_bf=41;
a_grid_bf=5*linspace(0,1,n_a_bf)'.^3;
pi_z_cpu=gather(pi_z);
z_grid_cpu=gather(z_grid);
d_grid_cpu=gather(d_grid);
U_j=zeros(n_d*n_a_bf,n_a_bf,n_z);
T_j=zeros(n_d*n_a_bf,n_a_bf,n_z);
V_bf=zeros(n_a_bf,n_z,N_j);
P_bf=zeros(n_a_bf,n_z,N_j); % joint (d,aprime) index, d fastest
for jj=N_j:-1:1
    for z_c=1:n_z
        for a_c=1:n_a_bf
            for ap_c=1:n_a_bf
                for d_c=1:n_d
                    U_j(d_c+n_d*(ap_c-1),a_c,z_c)=ReturnFn_z(d_grid_cpu(d_c),a_grid_bf(ap_c),a_grid_bf(a_c),z_grid_cpu(z_c),Params.r,Params.w,Params.kappa_j(jj),Params.sigma,Params.eta,Params.varphi,Params.agej(jj),Params.Jr,Params.pension);
                    T_j(d_c+n_d*(ap_c-1),a_c,z_c)=TemptationFn_z(d_grid_cpu(d_c),a_grid_bf(ap_c),a_grid_bf(a_c),z_grid_cpu(z_c),Params.lambdaGP,Params.shiftGP,Params.r,Params.w,Params.kappa_j(jj),Params.sigma,Params.agej(jj),Params.Jr,Params.pension);
                end
            end
        end
    end
    if jj==N_j
        EVbig=zeros(n_d*n_a_bf,1,n_z);
    else
        EVmat=V_bf(:,:,jj+1)*pi_z_cpu'; % E[V_{j+1}(a',z')|z], (aprime,z)
        EVbig=reshape(repelem(EVmat,n_d,1),[n_d*n_a_bf,1,n_z]); % expand to the joint (d,aprime) rows
    end
    MostT=max(T_j,[],1); % full-joint-choice-set most-tempting, (1,a,z)
    entireRHS=U_j+T_j+Params.beta*EVbig;
    [Vtmp,Ptmp]=max(entireRHS,[],1);
    V_bf(:,:,jj)=reshape(Vtmp-MostT,[n_a_bf,n_z]);
    P_bf(:,:,jj)=reshape(Ptmp,[n_a_bf,n_z]);
end
[Vtk1,Ptk1]=ValueFnIter_Case1_FHorz(n_d,n_a_bf,n_z,N_j,d_grid,a_grid_bf,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfgp1);
Ptk1joint=squeeze(Ptk1(1,:,:,:)+n_d*(Ptk1(2,:,:,:)-1));
fprintf('Cross test 4 (z_noe, plain): toolkit vs hand-rolled brute force, this should be zero: %.3e \n',max(abs(gather(Vtk1(:))-V_bf(:))))
fprintf('Cross test 4 (z_noe, plain): toolkit vs hand-rolled brute force, this should be zero: %.3e \n',max(abs(gather(Ptk1joint(:))-P_bf(:))))
[Vtk2,Ptk2]=ValueFnIter_Case1_FHorz(n_d,n_a_bf,n_z,N_j,d_grid,a_grid_bf,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfgp2);
Ptk2joint=squeeze(Ptk2(1,:,:,:)+n_d*(Ptk2(2,:,:,:)-1));
fprintf('Cross test 4 (z_noe, DC): toolkit vs hand-rolled brute force, this should be zero: %.3e \n',max(abs(gather(Vtk2(:))-V_bf(:))))
fprintf('Cross test 4 (z_noe, DC): toolkit vs hand-rolled brute force, this should be zero: %.3e \n',max(abs(gather(Ptk2joint(:))-P_bf(:))))
clear U_j T_j V_bf P_bf Vtk1 Vtk2 Ptk1 Ptk2 Ptk1joint Ptk2joint

%% ==================== Cross test 5: age-dependent lambdaGP via V_Jplus1 ====================
% Temptation in the first half of life only: lambdaGP_j=0.1 for j<=half, 0 after. The second
% half must equal the standard solve; the first half must equal a short (N_j=half) GP solve
% with constant lambdaGP and V_Jplus1 taken from the standard solve. At all four tiers, so the
% new raws' V_Jplus1 branches get runtime coverage (they are otherwise only reached in-loop).
Params5=Params;
Params5.lambdaGP=[0.1*ones(1,half),zeros(1,N_j-half)];
Njs=half;
Paramsjs=Params;
Paramsjs.agej=Params.agej(1:Njs);
Paramsjs.kappa_j=Params.kappa_j(1:Njs);
Paramsjs.lambdaGP=0.1;
% (mewj is age-dependent but is only used for the agent distribution, which is not computed here)
% plain tier
[V5a,P5a]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params5,DiscountFactorParamNames,[],vfgp1);
V5a_r=reshape(V5a,[],N_j);
P5a_r=reshape(P5a,[],N_j);
Vstd_r=reshape(Vstd1,[],N_j);
Pstd_r=reshape(Pstd1,[],N_j);
temp=V5a_r(:,half+1:N_j)-Vstd_r(:,half+1:N_j);
fprintf('Cross test 5 (z_noe, plain): second half (lambda=0) vs standard, this should be zero: %.3e \n',max(abs(temp(:))))
temp=P5a_r(:,half+1:N_j)-Pstd_r(:,half+1:N_j);
fprintf('Cross test 5 (z_noe, plain): second half (lambda=0) vs standard, this should be zero: %.3e \n',max(abs(temp(:))))
vfgp1c=vfgp1;
vfgp1c.V_Jplus1=Vstd1(:,:,half+1);
[V5c,P5c]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Paramsjs,DiscountFactorParamNames,[],vfgp1c);
temp=V5a_r(:,1:half)-reshape(V5c,[],Njs);
fprintf('Cross test 5 (z_noe, plain): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
temp=P5a_r(:,1:half)-reshape(P5c,[],Njs);
fprintf('Cross test 5 (z_noe, plain): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
% DC tier
[V5a,P5a]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params5,DiscountFactorParamNames,[],vfgp2);
V5a_r=reshape(V5a,[],N_j);
P5a_r=reshape(P5a,[],N_j);
Vstd_r=reshape(Vstd2,[],N_j);
Pstd_r=reshape(Pstd2,[],N_j);
temp=V5a_r(:,half+1:N_j)-Vstd_r(:,half+1:N_j);
fprintf('Cross test 5 (z_noe, DC): second half (lambda=0) vs standard, this should be zero: %.3e \n',max(abs(temp(:))))
temp=P5a_r(:,half+1:N_j)-Pstd_r(:,half+1:N_j);
fprintf('Cross test 5 (z_noe, DC): second half (lambda=0) vs standard, this should be zero: %.3e \n',max(abs(temp(:))))
vfgp2c=vfgp2;
vfgp2c.V_Jplus1=Vstd2(:,:,half+1);
[V5c,P5c]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Paramsjs,DiscountFactorParamNames,[],vfgp2c);
temp=V5a_r(:,1:half)-reshape(V5c,[],Njs);
fprintf('Cross test 5 (z_noe, DC): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
temp=P5a_r(:,1:half)-reshape(P5c,[],Njs);
fprintf('Cross test 5 (z_noe, DC): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
% GI tier
[V5a,P5a]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params5,DiscountFactorParamNames,[],vfgp3);
V5a_r=reshape(V5a,[],N_j);
P5a_r=reshape(P5a,[],N_j);
Vstd_r=reshape(Vstd3,[],N_j);
Pstd_r=reshape(Pstd3,[],N_j);
temp=V5a_r(:,half+1:N_j)-Vstd_r(:,half+1:N_j);
fprintf('Cross test 5 (z_noe, GI): second half (lambda=0) vs standard, this should be zero: %.3e \n',max(abs(temp(:))))
temp=P5a_r(:,half+1:N_j)-Pstd_r(:,half+1:N_j);
fprintf('Cross test 5 (z_noe, GI): second half (lambda=0) vs standard, this should be zero: %.3e \n',max(abs(temp(:))))
vfgp3c=vfgp3;
vfgp3c.V_Jplus1=Vstd3(:,:,half+1);
[V5c,P5c]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Paramsjs,DiscountFactorParamNames,[],vfgp3c);
temp=V5a_r(:,1:half)-reshape(V5c,[],Njs);
fprintf('Cross test 5 (z_noe, GI): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
temp=P5a_r(:,1:half)-reshape(P5c,[],Njs);
fprintf('Cross test 5 (z_noe, GI): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
% DC+GI tier
[V5a,P5a]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params5,DiscountFactorParamNames,[],vfgp4);
V5a_r=reshape(V5a,[],N_j);
P5a_r=reshape(P5a,[],N_j);
Vstd_r=reshape(Vstd4,[],N_j);
Pstd_r=reshape(Pstd4,[],N_j);
temp=V5a_r(:,half+1:N_j)-Vstd_r(:,half+1:N_j);
fprintf('Cross test 5 (z_noe, DC+GI): second half (lambda=0) vs standard, this should be zero: %.3e \n',max(abs(temp(:))))
temp=P5a_r(:,half+1:N_j)-Pstd_r(:,half+1:N_j);
fprintf('Cross test 5 (z_noe, DC+GI): second half (lambda=0) vs standard, this should be zero: %.3e \n',max(abs(temp(:))))
vfgp4c=vfgp4;
vfgp4c.V_Jplus1=Vstd4(:,:,half+1);
[V5c,P5c]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Paramsjs,DiscountFactorParamNames,[],vfgp4c);
temp=V5a_r(:,1:half)-reshape(V5c,[],Njs);
fprintf('Cross test 5 (z_noe, DC+GI): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
temp=P5a_r(:,1:half)-reshape(P5c,[],Njs);
fprintf('Cross test 5 (z_noe, DC+GI): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
clear V5a V5c P5a P5c V5a_r P5a_r Vstd_r Pstd_r

%% ==================== Cross test 6: temptation weakly lowers welfare ====================
% Self-control cost is nonnegative, so V_GP <= V_standard everywhere (lambdaGP=0.1 baseline)
temp=max(V1a(:)-Vstd1(:));
fprintf('Cross test 6 (z_noe, plain): max(V_GP - V_standard), this should be weakly negative: %.3e \n',temp)

%%
output=struct(); % Not currently used for anything. Maybe will do so later.

end
