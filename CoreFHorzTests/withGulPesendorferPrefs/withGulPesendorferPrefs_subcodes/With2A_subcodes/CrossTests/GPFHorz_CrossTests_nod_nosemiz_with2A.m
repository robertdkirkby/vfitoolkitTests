function output=GPFHorz_CrossTests_nod_nosemiz_with2A(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)
% The Gul-Pesendorfer with2A cross tests (without d): the same tests run on the
% two-endogenous-state model, pinning each GP 2A raw ({plain, DC2A, GI2A, DC2A+GI2A}) to its
% exponential donor via the lambdaGP=0 limit, plus the constant-temptation, shift-invariance,
% V_Jplus1 and sign tests. The hand-rolled brute-force tests of the nod donor are NOT
% repeated here: the main tier already pinned the temptation logic against the hand-rolls,
% and the 2A raws' new risk is indexing, which the lambda=0 and constant-v anchors catch.
% Test 1: lambdaGP=0 (no temptation) equals standard preferences. Run for all four shock
%         combos (z_noe, noz_e, z_e, noz_noe) at all four tiers (plain, DC2A, GI2A,
%         DC2A+GI2A), so every raw gets a limiting-case anchor.
% Test 2: constant temptation (lambdaGP=0, shiftGP=0.7) equals standard preferences: a
%         constant v cancels exactly against the most-tempting term.
% Test 3: shift invariance: v and v+0.7 give identical V and Policy (the shift cancels).
% Test 4: age-dependent lambdaGP (temptation in the first half of life only): second half
%         equals standard, first half equals a V_Jplus1 short solve. At all four tiers (this
%         also gives the new raws' V_Jplus1 branches runtime coverage).
% Test 5: sign check: temptation weakly lowers welfare, V_GP <= V_standard everywhere.

n_d=0;
d_grid=[];

ReturnFn_z=@(a1prime,a2prime,a1,a2,z,r,w,kappa_j,sigma,agej,Jr,pension,phi1,phi2) ReturnFn_nod_z_noe_nosemiz_with2A(a1prime,a2prime,a1,a2,z,r,w,kappa_j,sigma,agej,Jr,pension,phi1,phi2);
ReturnFn_e=@(a1prime,a2prime,a1,a2,e,r,w,kappa_j,sigma,agej,Jr,pension,phi1,phi2) ReturnFn_nod_noz_e_nosemiz_with2A(a1prime,a2prime,a1,a2,e,r,w,kappa_j,sigma,agej,Jr,pension,phi1,phi2);
ReturnFn_ze=@(a1prime,a2prime,a1,a2,z,e,r,w,kappa_j,sigma,agej,Jr,pension,phi1,phi2) ReturnFn_nod_z_e_nosemiz_with2A(a1prime,a2prime,a1,a2,z,e,r,w,kappa_j,sigma,agej,Jr,pension,phi1,phi2);
ReturnFn_none=@(a1prime,a2prime,a1,a2,r,w,kappa_j,sigma,agej,Jr,pension,phi1,phi2) ReturnFn_nod_noz_noe_nosemiz_with2A(a1prime,a2prime,a1,a2,r,w,kappa_j,sigma,agej,Jr,pension,phi1,phi2);

TemptationFn_z=@(a1prime,a2prime,a1,a2,z,lambdaGP,shiftGP,r,w,kappa_j,sigma,agej,Jr,pension) GPTemptationFn_nod_z_noe_nosemiz_with2A(a1prime,a2prime,a1,a2,z,lambdaGP,shiftGP,r,w,kappa_j,sigma,agej,Jr,pension);
TemptationFn_e=@(a1prime,a2prime,a1,a2,e,lambdaGP,shiftGP,r,w,kappa_j,sigma,agej,Jr,pension) GPTemptationFn_nod_noz_e_nosemiz_with2A(a1prime,a2prime,a1,a2,e,lambdaGP,shiftGP,r,w,kappa_j,sigma,agej,Jr,pension);
TemptationFn_ze=@(a1prime,a2prime,a1,a2,z,e,lambdaGP,shiftGP,r,w,kappa_j,sigma,agej,Jr,pension) GPTemptationFn_nod_z_e_nosemiz_with2A(a1prime,a2prime,a1,a2,z,e,lambdaGP,shiftGP,r,w,kappa_j,sigma,agej,Jr,pension);
TemptationFn_none=@(a1prime,a2prime,a1,a2,lambdaGP,shiftGP,r,w,kappa_j,sigma,agej,Jr,pension) GPTemptationFn_nod_noz_noe_nosemiz_with2A(a1prime,a2prime,a1,a2,lambdaGP,shiftGP,r,w,kappa_j,sigma,agej,Jr,pension);

n_e=vfoptionsbaseline.n_e;
e_grid=vfoptionsbaseline.e_grid;
pi_e=vfoptionsbaseline.pi_e;

half=N_j/2; % for cross test 4 (N_j=20)

% The four tiers as vfoptions add-ons (fields copied onto each solve's vfoptions below)
% plain: nothing; DC: divideandconquer=1; GI: gridinterplayer=1,ngridinterp=5; DC+GI: both

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
fprintf('Cross test 1 (z_noe with2A, plain): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(V1(:)-Vstd1(:))))
fprintf('Cross test 1 (z_noe with2A, plain): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(P1(:)-Pstd1(:))))
[V2,P2]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params0,DiscountFactorParamNames,[],vfgp2);
fprintf('Cross test 1 (z_noe with2A, DC): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(V2(:)-Vstd2(:))))
fprintf('Cross test 1 (z_noe with2A, DC): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(P2(:)-Pstd2(:))))
[V3,P3]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params0,DiscountFactorParamNames,[],vfgp3);
fprintf('Cross test 1 (z_noe with2A, GI): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(V3(:)-Vstd3(:))))
fprintf('Cross test 1 (z_noe with2A, GI): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(P3(:)-Pstd3(:))))
[V4,P4]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params0,DiscountFactorParamNames,[],vfgp4);
fprintf('Cross test 1 (z_noe with2A, DC+GI): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(V4(:)-Vstd4(:))))
fprintf('Cross test 1 (z_noe with2A, DC+GI): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(P4(:)-Pstd4(:))))
clear V1 V2 V3 V4 P1 P2 P3 P4
% (Vstd1..4/Pstd1..4 for the z_noe combo are kept: tests 2, 4 and 5 reuse them)

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
fprintf('Cross test 1 (noz_e with2A, plain): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(V1(:)-Vstd(:))))
fprintf('Cross test 1 (noz_e with2A, plain): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(P1(:)-Pstd(:))))
[Vstd,Pstd]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params,DiscountFactorParamNames,[],vfstd2e);
[V1,P1]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params0,DiscountFactorParamNames,[],vfgp2e);
fprintf('Cross test 1 (noz_e with2A, DC): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(V1(:)-Vstd(:))))
fprintf('Cross test 1 (noz_e with2A, DC): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(P1(:)-Pstd(:))))
[Vstd,Pstd]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params,DiscountFactorParamNames,[],vfstd3e);
[V1,P1]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params0,DiscountFactorParamNames,[],vfgp3e);
fprintf('Cross test 1 (noz_e with2A, GI): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(V1(:)-Vstd(:))))
fprintf('Cross test 1 (noz_e with2A, GI): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(P1(:)-Pstd(:))))
[Vstd,Pstd]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params,DiscountFactorParamNames,[],vfstd4e);
[V1,P1]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params0,DiscountFactorParamNames,[],vfgp4e);
fprintf('Cross test 1 (noz_e with2A, DC+GI): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(V1(:)-Vstd(:))))
fprintf('Cross test 1 (noz_e with2A, DC+GI): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(P1(:)-Pstd(:))))
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
fprintf('Cross test 1 (z_e with2A, plain): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(V1(:)-Vstd(:))))
fprintf('Cross test 1 (z_e with2A, plain): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(P1(:)-Pstd(:))))
[Vstd,Pstd]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_ze,Params,DiscountFactorParamNames,[],vfstd2ze);
[V1,P1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_ze,Params0,DiscountFactorParamNames,[],vfgp2ze);
fprintf('Cross test 1 (z_e with2A, DC): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(V1(:)-Vstd(:))))
fprintf('Cross test 1 (z_e with2A, DC): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(P1(:)-Pstd(:))))
[Vstd,Pstd]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_ze,Params,DiscountFactorParamNames,[],vfstd3ze);
[V1,P1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_ze,Params0,DiscountFactorParamNames,[],vfgp3ze);
fprintf('Cross test 1 (z_e with2A, GI): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(V1(:)-Vstd(:))))
fprintf('Cross test 1 (z_e with2A, GI): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(P1(:)-Pstd(:))))
[Vstd,Pstd]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_ze,Params,DiscountFactorParamNames,[],vfstd4ze);
[V1,P1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_ze,Params0,DiscountFactorParamNames,[],vfgp4ze);
fprintf('Cross test 1 (z_e with2A, DC+GI): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(V1(:)-Vstd(:))))
fprintf('Cross test 1 (z_e with2A, DC+GI): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(P1(:)-Pstd(:))))
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
fprintf('Cross test 1 (noz_noe with2A, plain): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(V1(:)-Vstd(:))))
fprintf('Cross test 1 (noz_noe with2A, plain): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(P1(:)-Pstd(:))))
[Vstd,Pstd]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_none,Params,DiscountFactorParamNames,[],vfstd2n);
[V1,P1]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_none,Params0,DiscountFactorParamNames,[],vfgp2n);
fprintf('Cross test 1 (noz_noe with2A, DC): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(V1(:)-Vstd(:))))
fprintf('Cross test 1 (noz_noe with2A, DC): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(P1(:)-Pstd(:))))
[Vstd,Pstd]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_none,Params,DiscountFactorParamNames,[],vfstd3n);
[V1,P1]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_none,Params0,DiscountFactorParamNames,[],vfgp3n);
fprintf('Cross test 1 (noz_noe with2A, GI): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(V1(:)-Vstd(:))))
fprintf('Cross test 1 (noz_noe with2A, GI): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(P1(:)-Pstd(:))))
[Vstd,Pstd]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_none,Params,DiscountFactorParamNames,[],vfstd4n);
[V1,P1]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_none,Params0,DiscountFactorParamNames,[],vfgp4n);
fprintf('Cross test 1 (noz_noe with2A, DC+GI): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(V1(:)-Vstd(:))))
fprintf('Cross test 1 (noz_noe with2A, DC+GI): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(P1(:)-Pstd(:))))
clear V1 P1 Vstd Pstd

%% ==================== Cross test 2: constant temptation equals standard ====================
% v = 0.7 on the feasible set (lambdaGP=0, shiftGP=0.7): the constant cancels exactly against
% the most-tempting term, so this equals standard preferences DESPITE the temptation machinery
% being fully engaged (nonzero MostTempting). Catches a dropped or double-counted subtraction.
Params2=Params;
Params2.lambdaGP=0;
Params2.shiftGP=0.7;
[V1,P1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params2,DiscountFactorParamNames,[],vfgp1);
fprintf('Cross test 2 (z_noe with2A, plain): constant temptation equals standard, this should be zero: %.3e \n',max(abs(V1(:)-Vstd1(:))))
fprintf('Cross test 2 (z_noe with2A, plain): constant temptation equals standard, this should be zero: %.3e \n',max(abs(P1(:)-Pstd1(:))))
[V2,P2]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params2,DiscountFactorParamNames,[],vfgp2);
fprintf('Cross test 2 (z_noe with2A, DC): constant temptation equals standard, this should be zero: %.3e \n',max(abs(V2(:)-Vstd2(:))))
fprintf('Cross test 2 (z_noe with2A, DC): constant temptation equals standard, this should be zero: %.3e \n',max(abs(P2(:)-Pstd2(:))))
[V3,P3]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params2,DiscountFactorParamNames,[],vfgp3);
fprintf('Cross test 2 (z_noe with2A, GI): constant temptation equals standard, this should be zero: %.3e \n',max(abs(V3(:)-Vstd3(:))))
fprintf('Cross test 2 (z_noe with2A, GI): constant temptation equals standard, this should be zero: %.3e \n',max(abs(P3(:)-Pstd3(:))))
[V4,P4]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params2,DiscountFactorParamNames,[],vfgp4);
fprintf('Cross test 2 (z_noe with2A, DC+GI): constant temptation equals standard, this should be zero: %.3e \n',max(abs(V4(:)-Vstd4(:))))
fprintf('Cross test 2 (z_noe with2A, DC+GI): constant temptation equals standard, this should be zero: %.3e \n',max(abs(P4(:)-Pstd4(:))))
clear V1 V2 V3 V4 P1 P2 P3 P4

%% ==================== Cross test 3: shift invariance ====================
% v and v+0.7 (lambdaGP=0.1, shiftGP 0 vs 0.7) give identical V and Policy
Params3a=Params; % lambdaGP=0.1, shiftGP=0 (the baseline)
Params3b=Params;
Params3b.shiftGP=0.7;
[V1a,P1a]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params3a,DiscountFactorParamNames,[],vfgp1);
[V1b,P1b]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params3b,DiscountFactorParamNames,[],vfgp1);
fprintf('Cross test 3 (z_noe with2A, plain): shift invariance, this should be zero: %.3e \n',max(abs(V1a(:)-V1b(:))))
fprintf('Cross test 3 (z_noe with2A, plain): shift invariance, this should be zero: %.3e \n',max(abs(P1a(:)-P1b(:))))
[V2a,P2a]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params3a,DiscountFactorParamNames,[],vfgp2);
[V2b,P2b]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params3b,DiscountFactorParamNames,[],vfgp2);
fprintf('Cross test 3 (z_noe with2A, DC): shift invariance, this should be zero: %.3e \n',max(abs(V2a(:)-V2b(:))))
fprintf('Cross test 3 (z_noe with2A, DC): shift invariance, this should be zero: %.3e \n',max(abs(P2a(:)-P2b(:))))
[V3a,P3a]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params3a,DiscountFactorParamNames,[],vfgp3);
[V3b,P3b]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params3b,DiscountFactorParamNames,[],vfgp3);
fprintf('Cross test 3 (z_noe with2A, GI): shift invariance, this should be zero: %.3e \n',max(abs(V3a(:)-V3b(:))))
fprintf('Cross test 3 (z_noe with2A, GI): shift invariance, this should be zero: %.3e \n',max(abs(P3a(:)-P3b(:))))
[V4a,P4a]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params3a,DiscountFactorParamNames,[],vfgp4);
[V4b,P4b]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params3b,DiscountFactorParamNames,[],vfgp4);
fprintf('Cross test 3 (z_noe with2A, DC+GI): shift invariance, this should be zero: %.3e \n',max(abs(V4a(:)-V4b(:))))
fprintf('Cross test 3 (z_noe with2A, DC+GI): shift invariance, this should be zero: %.3e \n',max(abs(P4a(:)-P4b(:))))
clear V2a V2b V3a V3b V4a V4b P2a P2b P3a P3b P4a P4b
% (V1a/P1a, the plain-tier GP solve at the baseline lambdaGP=0.1, is kept for test 5)

%% ==================== Cross test 4: age-dependent lambdaGP via V_Jplus1 ====================
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
fprintf('Cross test 4 (z_noe with2A, plain): second half (lambda=0) vs standard, this should be zero: %.3e \n',max(abs(temp(:))))
temp=P5a_r(:,half+1:N_j)-Pstd_r(:,half+1:N_j);
fprintf('Cross test 4 (z_noe with2A, plain): second half (lambda=0) vs standard, this should be zero: %.3e \n',max(abs(temp(:))))
vfgp1c=vfgp1;
vfgp1c.V_Jplus1=Vstd1(:,:,:,half+1);
[V5c,P5c]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Paramsjs,DiscountFactorParamNames,[],vfgp1c);
temp=V5a_r(:,1:half)-reshape(V5c,[],Njs);
fprintf('Cross test 4 (z_noe with2A, plain): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
temp=P5a_r(:,1:half)-reshape(P5c,[],Njs);
fprintf('Cross test 4 (z_noe with2A, plain): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
% DC tier
[V5a,P5a]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params5,DiscountFactorParamNames,[],vfgp2);
V5a_r=reshape(V5a,[],N_j);
P5a_r=reshape(P5a,[],N_j);
Vstd_r=reshape(Vstd2,[],N_j);
Pstd_r=reshape(Pstd2,[],N_j);
temp=V5a_r(:,half+1:N_j)-Vstd_r(:,half+1:N_j);
fprintf('Cross test 4 (z_noe with2A, DC): second half (lambda=0) vs standard, this should be zero: %.3e \n',max(abs(temp(:))))
temp=P5a_r(:,half+1:N_j)-Pstd_r(:,half+1:N_j);
fprintf('Cross test 4 (z_noe with2A, DC): second half (lambda=0) vs standard, this should be zero: %.3e \n',max(abs(temp(:))))
vfgp2c=vfgp2;
vfgp2c.V_Jplus1=Vstd2(:,:,:,half+1);
[V5c,P5c]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Paramsjs,DiscountFactorParamNames,[],vfgp2c);
temp=V5a_r(:,1:half)-reshape(V5c,[],Njs);
fprintf('Cross test 4 (z_noe with2A, DC): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
temp=P5a_r(:,1:half)-reshape(P5c,[],Njs);
fprintf('Cross test 4 (z_noe with2A, DC): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
% GI tier
[V5a,P5a]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params5,DiscountFactorParamNames,[],vfgp3);
V5a_r=reshape(V5a,[],N_j);
P5a_r=reshape(P5a,[],N_j);
Vstd_r=reshape(Vstd3,[],N_j);
Pstd_r=reshape(Pstd3,[],N_j);
temp=V5a_r(:,half+1:N_j)-Vstd_r(:,half+1:N_j);
fprintf('Cross test 4 (z_noe with2A, GI): second half (lambda=0) vs standard, this should be zero: %.3e \n',max(abs(temp(:))))
temp=P5a_r(:,half+1:N_j)-Pstd_r(:,half+1:N_j);
fprintf('Cross test 4 (z_noe with2A, GI): second half (lambda=0) vs standard, this should be zero: %.3e \n',max(abs(temp(:))))
vfgp3c=vfgp3;
vfgp3c.V_Jplus1=Vstd3(:,:,:,half+1);
[V5c,P5c]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Paramsjs,DiscountFactorParamNames,[],vfgp3c);
temp=V5a_r(:,1:half)-reshape(V5c,[],Njs);
fprintf('Cross test 4 (z_noe with2A, GI): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
temp=P5a_r(:,1:half)-reshape(P5c,[],Njs);
fprintf('Cross test 4 (z_noe with2A, GI): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
% DC+GI tier
[V5a,P5a]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params5,DiscountFactorParamNames,[],vfgp4);
V5a_r=reshape(V5a,[],N_j);
P5a_r=reshape(P5a,[],N_j);
Vstd_r=reshape(Vstd4,[],N_j);
Pstd_r=reshape(Pstd4,[],N_j);
temp=V5a_r(:,half+1:N_j)-Vstd_r(:,half+1:N_j);
fprintf('Cross test 4 (z_noe with2A, DC+GI): second half (lambda=0) vs standard, this should be zero: %.3e \n',max(abs(temp(:))))
temp=P5a_r(:,half+1:N_j)-Pstd_r(:,half+1:N_j);
fprintf('Cross test 4 (z_noe with2A, DC+GI): second half (lambda=0) vs standard, this should be zero: %.3e \n',max(abs(temp(:))))
vfgp4c=vfgp4;
vfgp4c.V_Jplus1=Vstd4(:,:,:,half+1);
[V5c,P5c]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Paramsjs,DiscountFactorParamNames,[],vfgp4c);
temp=V5a_r(:,1:half)-reshape(V5c,[],Njs);
fprintf('Cross test 4 (z_noe with2A, DC+GI): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
temp=P5a_r(:,1:half)-reshape(P5c,[],Njs);
fprintf('Cross test 4 (z_noe with2A, DC+GI): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
clear V5a V5c P5a P5c V5a_r P5a_r Vstd_r Pstd_r

%% ==================== Cross test 5: temptation weakly lowers welfare ====================
% Self-control cost is nonnegative, so V_GP <= V_standard everywhere (lambdaGP=0.1 baseline)
temp=max(V1a(:)-Vstd1(:));
fprintf('Cross test 5 (z_noe with2A, plain): max(V_GP - V_standard), this should be weakly negative: %.3e \n',temp)

%%
output=struct(); % Not currently used for anything. Maybe will do so later.

end
