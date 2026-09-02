function output=AmbRiskyAsset_CrossTests_d1(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,n_a_withA1,a_grid_withA1)
% RiskyAsset + AmbiguityAversion cross tests. u is treated as ambiguity, not risk: every solve declares
% ambiguity_pi_u. Every comparison is against exoticpreferences='None' — exact zeros, no tolerances.
% T1: identical priors = vNM (u&z on the noa1 tier; u&z at ALL FOUR methods on the withA1 tier; u&e noa1)
% T2: 3-pi-vs-9-pi duplicated u-priors (noz_noe model: pure return-distribution ambiguity)
% T3u a/b: an unambiguously worse pi_u binds, either slot — the ambiguity-averse investor behaves as if
%          facing the worst return distribution (exact: aprime increasing in u, V increasing in assets)
% T3z a/b: an unambiguously worse pi_z binds, either slot (u-priors replicated baseline)
% T4: age-varying n_ambiguity via V_Jplus1 (noa1 base; withA1 at all four methods)
% Plus: the pi_u-consistency warning stays silent/fires as appropriate.

ReturnFn_none=@(h,savings,a,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension) ReturnFn_d1_noz_noe_nosemiz(h,savings,a,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension);
ReturnFn_z=@(h,savings,a,z,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension) ReturnFn_d1_z_noe_nosemiz(h,savings,a,z,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension);
ReturnFn_e=@(h,savings,a,e,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension) ReturnFn_d1_noz_e_nosemiz(h,savings,a,e,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension);
ReturnFn_zA1=@(h,savings,a1prime,a1,a2,z,r,w,kappa_j,sigma,eta,varphi,r_a1,agej,Jr,pension) ReturnFn_d1_z_noe_nosemiz_withA1(h,savings,a1prime,a1,a2,z,r,w,kappa_j,sigma,eta,varphi,r_a1,agej,Jr,pension);

% The three baseline priors for each shock (prior 1 is the regular pi)
ambiguity_pi_u=vfoptionsbaseline.ambiguity_pi_u; % [n_u,3]
ambiguity_pi_z=vfoptionsbaseline.ambiguity_pi_z; % [n_z,n_z,3]
ambiguity_pi_e=vfoptionsbaseline.ambiguity_pi_e; % [n_e,3]
n_e=vfoptionsbaseline.n_e;
e_grid=vfoptionsbaseline.e_grid;
pi_e=vfoptionsbaseline.pi_e;
n_u=vfoptionsbaseline.n_u;
pi_u=vfoptionsbaseline.pi_u;

% riskyasset vfoptions shared by every solve in this file
rvfoptions=struct();
rvfoptions.riskyasset=1;
rvfoptions.refine_d=[1,1,1];
rvfoptions.aprimeFn=vfoptionsbaseline.aprimeFn;
rvfoptions.n_u=vfoptionsbaseline.n_u;
rvfoptions.u_grid=vfoptionsbaseline.u_grid;
rvfoptions.pi_u=vfoptionsbaseline.pi_u;

half=N_j/2; % for cross test 4 (N_j=20)

%% Cross test 1: three identical priors (every declared shock) = vNM
% (a) u & z ambiguous, noa1 tier, base method
[Vstdz,Policystdz]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],rvfoptions);
vfoptions_t1=rvfoptions;
vfoptions_t1.exoticpreferences='AmbiguityAversion';
vfoptions_t1.n_ambiguity=3;
vfoptions_t1.ambiguity_pi_u=repmat(pi_u,[1,3]);
vfoptions_t1.ambiguity_pi_z=cat(3,pi_z,pi_z,pi_z);
[V1,Policy1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfoptions_t1);
fprintf('RiskyAsset cross test 1 (u&z, noa1): three identical priors, this should be zero: %.3e \n',max(abs(V1(:)-Vstdz(:))))
fprintf('RiskyAsset cross test 1 (u&z, noa1): three identical priors, this should be zero: %.3e \n',max(abs(Policy1(:)-Policystdz(:))))
clear V1 Policy1
% (b) u & z ambiguous, withA1 tier, at all four methods
vfoptions_none1=rvfoptions;
[VstdA1,PolicystdA1]=ValueFnIter_Case1_FHorz(n_d,n_a_withA1,n_z,N_j,d_grid,a_grid_withA1,z_grid,pi_z,ReturnFn_zA1,Params,DiscountFactorParamNames,[],vfoptions_none1);
vfoptions_t11=vfoptions_none1;
vfoptions_t11.exoticpreferences='AmbiguityAversion';
vfoptions_t11.n_ambiguity=3;
vfoptions_t11.ambiguity_pi_u=repmat(pi_u,[1,3]);
vfoptions_t11.ambiguity_pi_z=cat(3,pi_z,pi_z,pi_z);
[V11,Policy11]=ValueFnIter_Case1_FHorz(n_d,n_a_withA1,n_z,N_j,d_grid,a_grid_withA1,z_grid,pi_z,ReturnFn_zA1,Params,DiscountFactorParamNames,[],vfoptions_t11);
fprintf('RiskyAsset cross test 1 (u&z, withA1, base): three identical priors, this should be zero: %.3e \n',max(abs(V11(:)-VstdA1(:))))
fprintf('RiskyAsset cross test 1 (u&z, withA1, base): three identical priors, this should be zero: %.3e \n',max(abs(Policy11(:)-PolicystdA1(:))))
clear V11 Policy11
vfoptions_none2=rvfoptions;
vfoptions_none2.divideandconquer=1;
[VstdA2,PolicystdA2]=ValueFnIter_Case1_FHorz(n_d,n_a_withA1,n_z,N_j,d_grid,a_grid_withA1,z_grid,pi_z,ReturnFn_zA1,Params,DiscountFactorParamNames,[],vfoptions_none2);
vfoptions_t12=vfoptions_none2;
vfoptions_t12.exoticpreferences='AmbiguityAversion';
vfoptions_t12.n_ambiguity=3;
vfoptions_t12.ambiguity_pi_u=repmat(pi_u,[1,3]);
vfoptions_t12.ambiguity_pi_z=cat(3,pi_z,pi_z,pi_z);
[V12,Policy12]=ValueFnIter_Case1_FHorz(n_d,n_a_withA1,n_z,N_j,d_grid,a_grid_withA1,z_grid,pi_z,ReturnFn_zA1,Params,DiscountFactorParamNames,[],vfoptions_t12);
fprintf('RiskyAsset cross test 1 (u&z, withA1, DC): three identical priors, this should be zero: %.3e \n',max(abs(V12(:)-VstdA2(:))))
fprintf('RiskyAsset cross test 1 (u&z, withA1, DC): three identical priors, this should be zero: %.3e \n',max(abs(Policy12(:)-PolicystdA2(:))))
clear V12 Policy12
vfoptions_none3=rvfoptions;
vfoptions_none3.gridinterplayer=1;
vfoptions_none3.ngridinterp=5;
[VstdA3,PolicystdA3]=ValueFnIter_Case1_FHorz(n_d,n_a_withA1,n_z,N_j,d_grid,a_grid_withA1,z_grid,pi_z,ReturnFn_zA1,Params,DiscountFactorParamNames,[],vfoptions_none3);
vfoptions_t13=vfoptions_none3;
vfoptions_t13.exoticpreferences='AmbiguityAversion';
vfoptions_t13.n_ambiguity=3;
vfoptions_t13.ambiguity_pi_u=repmat(pi_u,[1,3]);
vfoptions_t13.ambiguity_pi_z=cat(3,pi_z,pi_z,pi_z);
[V13,Policy13]=ValueFnIter_Case1_FHorz(n_d,n_a_withA1,n_z,N_j,d_grid,a_grid_withA1,z_grid,pi_z,ReturnFn_zA1,Params,DiscountFactorParamNames,[],vfoptions_t13);
fprintf('RiskyAsset cross test 1 (u&z, withA1, GI): three identical priors, this should be zero: %.3e \n',max(abs(V13(:)-VstdA3(:))))
fprintf('RiskyAsset cross test 1 (u&z, withA1, GI): three identical priors, this should be zero: %.3e \n',max(abs(Policy13(:)-PolicystdA3(:))))
clear V13 Policy13
vfoptions_none4=rvfoptions;
vfoptions_none4.divideandconquer=1;
vfoptions_none4.gridinterplayer=1;
vfoptions_none4.ngridinterp=5;
[VstdA4,PolicystdA4]=ValueFnIter_Case1_FHorz(n_d,n_a_withA1,n_z,N_j,d_grid,a_grid_withA1,z_grid,pi_z,ReturnFn_zA1,Params,DiscountFactorParamNames,[],vfoptions_none4);
vfoptions_t14=vfoptions_none4;
vfoptions_t14.exoticpreferences='AmbiguityAversion';
vfoptions_t14.n_ambiguity=3;
vfoptions_t14.ambiguity_pi_u=repmat(pi_u,[1,3]);
vfoptions_t14.ambiguity_pi_z=cat(3,pi_z,pi_z,pi_z);
[V14,Policy14]=ValueFnIter_Case1_FHorz(n_d,n_a_withA1,n_z,N_j,d_grid,a_grid_withA1,z_grid,pi_z,ReturnFn_zA1,Params,DiscountFactorParamNames,[],vfoptions_t14);
fprintf('RiskyAsset cross test 1 (u&z, withA1, DC+GI): three identical priors, this should be zero: %.3e \n',max(abs(V14(:)-VstdA4(:))))
fprintf('RiskyAsset cross test 1 (u&z, withA1, DC+GI): three identical priors, this should be zero: %.3e \n',max(abs(Policy14(:)-PolicystdA4(:))))
clear V14 Policy14
% (c) u & e ambiguous, noa1 tier, base method
vfoptions_nonee=rvfoptions;
vfoptions_nonee.n_e=n_e;
vfoptions_nonee.e_grid=e_grid;
vfoptions_nonee.pi_e=pi_e;
[Vstde,Policystde]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params,DiscountFactorParamNames,[],vfoptions_nonee);
vfoptions_t1e=vfoptions_nonee;
vfoptions_t1e.exoticpreferences='AmbiguityAversion';
vfoptions_t1e.n_ambiguity=3;
vfoptions_t1e.ambiguity_pi_u=repmat(pi_u,[1,3]);
vfoptions_t1e.ambiguity_pi_e=[pi_e,pi_e,pi_e];
[V1e,Policy1e]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params,DiscountFactorParamNames,[],vfoptions_t1e);
fprintf('RiskyAsset cross test 1 (u&e, noa1): three identical priors, this should be zero: %.3e \n',max(abs(V1e(:)-Vstde(:))))
fprintf('RiskyAsset cross test 1 (u&e, noa1): three identical priors, this should be zero: %.3e \n',max(abs(Policy1e(:)-Policystde(:))))
clear V1e Policy1e Vstde Policystde

%% Cross test 2: duplicated u-priors, 3 pi vs 9 pi (pure return-distribution ambiguity, noz_noe noa1)
vfoptions_t2=rvfoptions;
vfoptions_t2.exoticpreferences='AmbiguityAversion';
vfoptions_t2.n_ambiguity=3;
vfoptions_t2.ambiguity_pi_u=ambiguity_pi_u;
[V2a,Policy2a]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_none,Params,DiscountFactorParamNames,[],vfoptions_t2);
vfoptions_t2b=rvfoptions;
vfoptions_t2b.exoticpreferences='AmbiguityAversion';
vfoptions_t2b.n_ambiguity=9;
vfoptions_t2b.ambiguity_pi_u=ambiguity_pi_u(:,[1,2,3,1,2,3,1,1,1]);
[V2b,Policy2b]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_none,Params,DiscountFactorParamNames,[],vfoptions_t2b);
fprintf('RiskyAsset cross test 2 (u): 3 pi vs 9 pi, this should be zero: %.3e \n',max(abs(V2a(:)-V2b(:))))
fprintf('RiskyAsset cross test 2 (u): 3 pi vs 9 pi, this should be zero: %.3e \n',max(abs(Policy2a(:)-Policy2b(:))))
clear V2a V2b Policy2a Policy2b

%% Cross tests 3u a/b: an unambiguously worse pi_u binds, in either prior slot (noz_noe noa1)
% The worse pi_u shifts half of the best (highest-return) u's probability onto the worst u
shiftu=0.5*pi_u(n_u);
pi_u_worse=pi_u;
pi_u_worse(n_u)=pi_u(n_u)-shiftu;
pi_u_worse(1)=pi_u(1)+shiftu;
% Standard preferences under the worse pi_u
vfoptions_nonew=rvfoptions;
vfoptions_nonew.pi_u=pi_u_worse;
[Vstdw,Policystdw]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_none,Params,DiscountFactorParamNames,[],vfoptions_nonew);
% 3a: worse pi_u in prior slot 1
vfoptions_t3=rvfoptions;
vfoptions_t3.exoticpreferences='AmbiguityAversion';
vfoptions_t3.n_ambiguity=2;
vfoptions_t3.ambiguity_pi_u=[pi_u_worse,pi_u];
[V3a,Policy3a]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_none,Params,DiscountFactorParamNames,[],vfoptions_t3);
fprintf('RiskyAsset cross test 3a (u): worse pi_u (slot 1) binds, this should be zero: %.3e \n',max(abs(V3a(:)-Vstdw(:))))
fprintf('RiskyAsset cross test 3a (u): worse pi_u (slot 1) binds, this should be zero: %.3e \n',max(abs(Policy3a(:)-Policystdw(:))))
% 3b: worse pi_u in prior slot 2 (ordering of the priors is irrelevant)
vfoptions_t3.ambiguity_pi_u=[pi_u,pi_u_worse];
[V3b,Policy3b]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_none,Params,DiscountFactorParamNames,[],vfoptions_t3);
fprintf('RiskyAsset cross test 3b (u): worse pi_u (slot 2) binds, this should be zero: %.3e \n',max(abs(V3b(:)-Vstdw(:))))
fprintf('RiskyAsset cross test 3b (u): worse pi_u (slot 2) binds, this should be zero: %.3e \n',max(abs(Policy3b(:)-Policystdw(:))))
clear V3a V3b Policy3a Policy3b Vstdw Policystdw

%% Cross tests 3z a/b: an unambiguously worse pi_z binds, in either prior slot (z noa1; u-priors replicated baseline)
shiftz=0.5*pi_z(:,n_z);
pi_z_worse=pi_z;
pi_z_worse(:,n_z)=pi_z(:,n_z)-shiftz;
pi_z_worse(:,1)=pi_z(:,1)+shiftz;
% Standard preferences under the worse pi_z
[Vstdwz,Policystdwz]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z_worse,ReturnFn_z,Params,DiscountFactorParamNames,[],rvfoptions);
% 3a: worse pi_z in prior slot 1
vfoptions_t3z=rvfoptions;
vfoptions_t3z.exoticpreferences='AmbiguityAversion';
vfoptions_t3z.n_ambiguity=2;
vfoptions_t3z.ambiguity_pi_u=repmat(pi_u,[1,2]);
vfoptions_t3z.ambiguity_pi_z=cat(3,pi_z_worse,pi_z);
[V3za,Policy3za]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfoptions_t3z);
fprintf('RiskyAsset cross test 3a (z): worse pi_z (slot 1) binds, this should be zero: %.3e \n',max(abs(V3za(:)-Vstdwz(:))))
fprintf('RiskyAsset cross test 3a (z): worse pi_z (slot 1) binds, this should be zero: %.3e \n',max(abs(Policy3za(:)-Policystdwz(:))))
% 3b: worse pi_z in prior slot 2 (ordering of the priors is irrelevant)
vfoptions_t3z.ambiguity_pi_z=cat(3,pi_z,pi_z_worse);
[V3zb,Policy3zb]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfoptions_t3z);
fprintf('RiskyAsset cross test 3b (z): worse pi_z (slot 2) binds, this should be zero: %.3e \n',max(abs(V3zb(:)-Vstdwz(:))))
fprintf('RiskyAsset cross test 3b (z): worse pi_z (slot 2) binds, this should be zero: %.3e \n',max(abs(Policy3zb(:)-Policystdwz(:))))
clear V3za V3zb Policy3za Policy3zb Vstdwz Policystdwz

%% Cross test 4: age-varying n_ambiguity via V_Jplus1 (three u- and z-priors in the first half, single
% baseline prior in the second half; second half is then vNM, and the first half must equal a short
% solve seeded with the vNM V at half+1)
% Slice jj of pi_z_J is the transition from period jj to jj+1, so the three z-priors sit in slices 1..half;
% ambiguity_pi_u is age-independent, and with n_ambiguity(jj)=1 in the second half only its first (baseline)
% column is read there.
ambiguity_pi_z_J=repmat(reshape(pi_z,[n_z,n_z,1,1]),[1,1,N_j,3]);
ambiguity_pi_z_J(:,:,1:half,2)=repmat(ambiguity_pi_z(:,:,2),[1,1,half]);
ambiguity_pi_z_J(:,:,1:half,3)=repmat(ambiguity_pi_z(:,:,3),[1,1,half]);
Njs=half;
Paramsjs=Params;
Paramsjs.agej=Params.agej(1:Njs);
Paramsjs.kappa_j=Params.kappa_j(1:Njs);
% (mewj is age-dependent but is only used for the agent distribution, which is not computed here)
% --- noa1 tier, base method ---
vfoptions_t4a=rvfoptions;
vfoptions_t4a.exoticpreferences='AmbiguityAversion';
vfoptions_t4a.n_ambiguity=[3*ones(1,half),ones(1,N_j-half)];
vfoptions_t4a.ambiguity_pi_u=ambiguity_pi_u;
vfoptions_t4a.ambiguity_pi_z_J=ambiguity_pi_z_J;
[V4a,Policy4a]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfoptions_t4a);
V4a_r=reshape(V4a,[],N_j);
P4a_r=reshape(Policy4a,[],N_j);
Vstd_r=reshape(Vstdz,[],N_j);
Pstd_r=reshape(Policystdz,[],N_j);
temp=V4a_r(:,half+1:N_j)-Vstd_r(:,half+1:N_j);
fprintf('RiskyAsset cross test 4 (noa1, base): second half (single prior) vs standard, this should be zero: %.3e \n',max(abs(temp(:))))
temp=P4a_r(:,half+1:N_j)-Pstd_r(:,half+1:N_j);
fprintf('RiskyAsset cross test 4 (noa1, base): second half (single prior) vs standard, this should be zero: %.3e \n',max(abs(temp(:))))
vfoptions_t4c=rvfoptions;
vfoptions_t4c.exoticpreferences='AmbiguityAversion';
vfoptions_t4c.n_ambiguity=3;
vfoptions_t4c.ambiguity_pi_u=ambiguity_pi_u;
vfoptions_t4c.ambiguity_pi_z=ambiguity_pi_z;
vfoptions_t4c.V_Jplus1=Vstdz(:,:,half+1);
[V4c,Policy4c]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Paramsjs,DiscountFactorParamNames,[],vfoptions_t4c);
temp=V4a_r(:,1:half)-reshape(V4c,[],Njs);
fprintf('RiskyAsset cross test 4 (noa1, base): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
temp=P4a_r(:,1:half)-reshape(Policy4c,[],Njs);
fprintf('RiskyAsset cross test 4 (noa1, base): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
clear V4a V4c Policy4a Policy4c Vstdz Policystdz
% --- withA1 tier, base method ---
vfoptions_t4a1=vfoptions_none1;
vfoptions_t4a1.exoticpreferences='AmbiguityAversion';
vfoptions_t4a1.n_ambiguity=[3*ones(1,half),ones(1,N_j-half)];
vfoptions_t4a1.ambiguity_pi_u=ambiguity_pi_u;
vfoptions_t4a1.ambiguity_pi_z_J=ambiguity_pi_z_J;
[V4a1,Policy4a1]=ValueFnIter_Case1_FHorz(n_d,n_a_withA1,n_z,N_j,d_grid,a_grid_withA1,z_grid,pi_z,ReturnFn_zA1,Params,DiscountFactorParamNames,[],vfoptions_t4a1);
V4a_r=reshape(V4a1,[],N_j);
P4a_r=reshape(Policy4a1,[],N_j);
Vstd_r=reshape(VstdA1,[],N_j);
Pstd_r=reshape(PolicystdA1,[],N_j);
temp=V4a_r(:,half+1:N_j)-Vstd_r(:,half+1:N_j);
fprintf('RiskyAsset cross test 4 (withA1, base): second half (single prior) vs standard, this should be zero: %.3e \n',max(abs(temp(:))))
temp=P4a_r(:,half+1:N_j)-Pstd_r(:,half+1:N_j);
fprintf('RiskyAsset cross test 4 (withA1, base): second half (single prior) vs standard, this should be zero: %.3e \n',max(abs(temp(:))))
vfoptions_t4c1=vfoptions_none1;
vfoptions_t4c1.exoticpreferences='AmbiguityAversion';
vfoptions_t4c1.n_ambiguity=3;
vfoptions_t4c1.ambiguity_pi_u=ambiguity_pi_u;
vfoptions_t4c1.ambiguity_pi_z=ambiguity_pi_z;
vfoptions_t4c1.V_Jplus1=VstdA1(:,:,:,half+1);
[V4c1,Policy4c1]=ValueFnIter_Case1_FHorz(n_d,n_a_withA1,n_z,Njs,d_grid,a_grid_withA1,z_grid,pi_z,ReturnFn_zA1,Paramsjs,DiscountFactorParamNames,[],vfoptions_t4c1);
temp=V4a_r(:,1:half)-reshape(V4c1,[],Njs);
fprintf('RiskyAsset cross test 4 (withA1, base): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
temp=P4a_r(:,1:half)-reshape(Policy4c1,[],Njs);
fprintf('RiskyAsset cross test 4 (withA1, base): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
clear V4a1 V4c1 Policy4a1 Policy4c1 VstdA1 PolicystdA1
% --- withA1 tier, DC method ---
vfoptions_t4a2=vfoptions_none2;
vfoptions_t4a2.exoticpreferences='AmbiguityAversion';
vfoptions_t4a2.n_ambiguity=[3*ones(1,half),ones(1,N_j-half)];
vfoptions_t4a2.ambiguity_pi_u=ambiguity_pi_u;
vfoptions_t4a2.ambiguity_pi_z_J=ambiguity_pi_z_J;
[V4a2,Policy4a2]=ValueFnIter_Case1_FHorz(n_d,n_a_withA1,n_z,N_j,d_grid,a_grid_withA1,z_grid,pi_z,ReturnFn_zA1,Params,DiscountFactorParamNames,[],vfoptions_t4a2);
V4a_r=reshape(V4a2,[],N_j);
P4a_r=reshape(Policy4a2,[],N_j);
Vstd_r=reshape(VstdA2,[],N_j);
Pstd_r=reshape(PolicystdA2,[],N_j);
temp=V4a_r(:,half+1:N_j)-Vstd_r(:,half+1:N_j);
fprintf('RiskyAsset cross test 4 (withA1, DC): second half (single prior) vs standard, this should be zero: %.3e \n',max(abs(temp(:))))
temp=P4a_r(:,half+1:N_j)-Pstd_r(:,half+1:N_j);
fprintf('RiskyAsset cross test 4 (withA1, DC): second half (single prior) vs standard, this should be zero: %.3e \n',max(abs(temp(:))))
vfoptions_t4c2=vfoptions_none2;
vfoptions_t4c2.exoticpreferences='AmbiguityAversion';
vfoptions_t4c2.n_ambiguity=3;
vfoptions_t4c2.ambiguity_pi_u=ambiguity_pi_u;
vfoptions_t4c2.ambiguity_pi_z=ambiguity_pi_z;
vfoptions_t4c2.V_Jplus1=VstdA2(:,:,:,half+1);
[V4c2,Policy4c2]=ValueFnIter_Case1_FHorz(n_d,n_a_withA1,n_z,Njs,d_grid,a_grid_withA1,z_grid,pi_z,ReturnFn_zA1,Paramsjs,DiscountFactorParamNames,[],vfoptions_t4c2);
temp=V4a_r(:,1:half)-reshape(V4c2,[],Njs);
fprintf('RiskyAsset cross test 4 (withA1, DC): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
temp=P4a_r(:,1:half)-reshape(Policy4c2,[],Njs);
fprintf('RiskyAsset cross test 4 (withA1, DC): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
clear V4a2 V4c2 Policy4a2 Policy4c2 VstdA2 PolicystdA2
% --- withA1 tier, GI method ---
vfoptions_t4a3=vfoptions_none3;
vfoptions_t4a3.exoticpreferences='AmbiguityAversion';
vfoptions_t4a3.n_ambiguity=[3*ones(1,half),ones(1,N_j-half)];
vfoptions_t4a3.ambiguity_pi_u=ambiguity_pi_u;
vfoptions_t4a3.ambiguity_pi_z_J=ambiguity_pi_z_J;
[V4a3,Policy4a3]=ValueFnIter_Case1_FHorz(n_d,n_a_withA1,n_z,N_j,d_grid,a_grid_withA1,z_grid,pi_z,ReturnFn_zA1,Params,DiscountFactorParamNames,[],vfoptions_t4a3);
V4a_r=reshape(V4a3,[],N_j);
P4a_r=reshape(Policy4a3,[],N_j);
Vstd_r=reshape(VstdA3,[],N_j);
Pstd_r=reshape(PolicystdA3,[],N_j);
temp=V4a_r(:,half+1:N_j)-Vstd_r(:,half+1:N_j);
fprintf('RiskyAsset cross test 4 (withA1, GI): second half (single prior) vs standard, this should be zero: %.3e \n',max(abs(temp(:))))
temp=P4a_r(:,half+1:N_j)-Pstd_r(:,half+1:N_j);
fprintf('RiskyAsset cross test 4 (withA1, GI): second half (single prior) vs standard, this should be zero: %.3e \n',max(abs(temp(:))))
vfoptions_t4c3=vfoptions_none3;
vfoptions_t4c3.exoticpreferences='AmbiguityAversion';
vfoptions_t4c3.n_ambiguity=3;
vfoptions_t4c3.ambiguity_pi_u=ambiguity_pi_u;
vfoptions_t4c3.ambiguity_pi_z=ambiguity_pi_z;
vfoptions_t4c3.V_Jplus1=VstdA3(:,:,:,half+1);
[V4c3,Policy4c3]=ValueFnIter_Case1_FHorz(n_d,n_a_withA1,n_z,Njs,d_grid,a_grid_withA1,z_grid,pi_z,ReturnFn_zA1,Paramsjs,DiscountFactorParamNames,[],vfoptions_t4c3);
temp=V4a_r(:,1:half)-reshape(V4c3,[],Njs);
fprintf('RiskyAsset cross test 4 (withA1, GI): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
temp=P4a_r(:,1:half)-reshape(Policy4c3,[],Njs);
fprintf('RiskyAsset cross test 4 (withA1, GI): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
clear V4a3 V4c3 Policy4a3 Policy4c3 VstdA3 PolicystdA3
% --- withA1 tier, DC+GI method ---
vfoptions_t4a4=vfoptions_none4;
vfoptions_t4a4.exoticpreferences='AmbiguityAversion';
vfoptions_t4a4.n_ambiguity=[3*ones(1,half),ones(1,N_j-half)];
vfoptions_t4a4.ambiguity_pi_u=ambiguity_pi_u;
vfoptions_t4a4.ambiguity_pi_z_J=ambiguity_pi_z_J;
[V4a4,Policy4a4]=ValueFnIter_Case1_FHorz(n_d,n_a_withA1,n_z,N_j,d_grid,a_grid_withA1,z_grid,pi_z,ReturnFn_zA1,Params,DiscountFactorParamNames,[],vfoptions_t4a4);
V4a_r=reshape(V4a4,[],N_j);
P4a_r=reshape(Policy4a4,[],N_j);
Vstd_r=reshape(VstdA4,[],N_j);
Pstd_r=reshape(PolicystdA4,[],N_j);
temp=V4a_r(:,half+1:N_j)-Vstd_r(:,half+1:N_j);
fprintf('RiskyAsset cross test 4 (withA1, DC+GI): second half (single prior) vs standard, this should be zero: %.3e \n',max(abs(temp(:))))
temp=P4a_r(:,half+1:N_j)-Pstd_r(:,half+1:N_j);
fprintf('RiskyAsset cross test 4 (withA1, DC+GI): second half (single prior) vs standard, this should be zero: %.3e \n',max(abs(temp(:))))
vfoptions_t4c4=vfoptions_none4;
vfoptions_t4c4.exoticpreferences='AmbiguityAversion';
vfoptions_t4c4.n_ambiguity=3;
vfoptions_t4c4.ambiguity_pi_u=ambiguity_pi_u;
vfoptions_t4c4.ambiguity_pi_z=ambiguity_pi_z;
vfoptions_t4c4.V_Jplus1=VstdA4(:,:,:,half+1);
[V4c4,Policy4c4]=ValueFnIter_Case1_FHorz(n_d,n_a_withA1,n_z,Njs,d_grid,a_grid_withA1,z_grid,pi_z,ReturnFn_zA1,Paramsjs,DiscountFactorParamNames,[],vfoptions_t4c4);
temp=V4a_r(:,1:half)-reshape(V4c4,[],Njs);
fprintf('RiskyAsset cross test 4 (withA1, DC+GI): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
temp=P4a_r(:,1:half)-reshape(Policy4c4,[],Njs);
fprintf('RiskyAsset cross test 4 (withA1, DC+GI): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
clear V4a4 V4c4 Policy4a4 Policy4c4 VstdA4 PolicystdA4

%% pi_u-consistency warning: the regular pi_u should be one of the u-priors
lastwarn(''); % reset
[Vw,Policyw]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_none,Params,DiscountFactorParamNames,[],vfoptions_t2);
warnmsg=lastwarn();
if isempty(warnmsg)
    fprintf('pi_u-consistency: no warning when pi_u is one of the priors :) \n')
else
    fprintf('pi_u-consistency: FAIL, got an unexpected warning: %s \n',warnmsg)
end
lastwarn(''); % reset
vfoptions_t2off=vfoptions_t2;
vfoptions_t2off.pi_u=0.5*pi_u+0.5*ones(n_u,1)/n_u; % not equal to any of the u-priors
[Vw,Policyw]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_none,Params,DiscountFactorParamNames,[],vfoptions_t2off);
warnmsg=lastwarn();
if isempty(warnmsg)
    fprintf('pi_u-consistency: FAIL, expected a warning (pi_u is not one of the priors) and got none \n')
else
    fprintf('pi_u-consistency: warning fires when pi_u is not one of the priors :) \n')
end

%%
output=struct(); % Not currently used for anything. Maybe will do so later.

end
