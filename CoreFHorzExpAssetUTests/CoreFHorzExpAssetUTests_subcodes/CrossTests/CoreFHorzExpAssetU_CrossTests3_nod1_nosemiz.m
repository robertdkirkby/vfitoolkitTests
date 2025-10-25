function output=CoreFHorzExpAssetU_CrossTests3_nod1_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)

% Just solve without z and without e
ReturnFn_twoendo=@(a1prime,a2prime,a1,a2,r,w,kappa_j,sigma,agej,Jr,pension)...
    ReturnFn_TwoEndo_nod1_noz_noe_nosemiz(a1prime,a2prime,a1,a2,r,w,kappa_j,sigma,agej,Jr,pension);

ReturnFn_none=@(d2,a1prime,a1,a2,r,w,kappa_j,sigma,agej,Jr,pension)...
    ReturnFn_nod1_noz_noe_nosemiz(d2,a1prime,a1,a2,r,w,kappa_j,sigma,agej,Jr,pension);

% Note: d2 and a2prime should end up just the same thing. The two endo
% problem is kind of stupid, but that it okay as it is just to check the
% compute.

aprimeFn=@(d2,a2,u) d2; % d2 is just a2prime
d_grid=a_grid(n_a(1)+1:end); % set d2 to a2prime (a2prime is just a2)
n_d=n_a(2); % set d2 to a2prime (a2prime is just a2)

%% Solving with just a single points for z with value 1 and prob 1 gives us same as no shocks
jequaloneDist_none=zeros([n_a],'gpuArray');
jequaloneDist_none(1,1)=1; % no assets

vfoptionsA.divideandconquer=1;
simoptionsA=struct();
[V0,Policy0]=ValueFnIter_Case1_FHorz(0,n_a,0,N_j,[],a_grid,[],[],ReturnFn_twoendo,Params,DiscountFactorParamNames,[],vfoptionsA);
StationaryDist0=StationaryDist_FHorz_Case1(jequaloneDist_none,AgeWeightParamNames,Policy0,0,n_a,0,N_j,[],Params,simoptionsA);

Policy0alt=[Policy0(2,:,:,:); Policy0(1,:,:,:)]; % swap order


% Experience asset u
vfoptionsB.experienceassetu=1;
simoptionsB.experienceassetu=1;
vfoptionsB.aprimeFn=aprimeFn;
simoptionsB.aprimeFn=aprimeFn;
simoptionsB.d_grid=d_grid;
simoptionsB.a_grid=a_grid;
vfoptionsB.n_u=1;
vfoptionsB.u_grid=1;
vfoptionsB.pi_u=1;
simoptionsB.n_u=1;
simoptionsB.u_grid=1;
simoptionsB.pi_u=1;
[V1,Policy1]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_none,Params,DiscountFactorParamNames,[],vfoptionsB);
StationaryDist1=StationaryDist_FHorz_Case1(jequaloneDist_none,AgeWeightParamNames,Policy1,n_d,n_a,0,N_j,1,Params,simoptionsB);

fprintf('Cross test 3: expasset is just a standard endo state, this should be zero: %2.8f \n',max(abs(V0(:)-V1(:))))
fprintf('Cross test 3: expasset is just a standard endo state, this should be zero: %2.8f \n',max(abs(Policy0alt(:)-Policy1(:))))
fprintf('Cross test 3: expasset is just a standard endo state, this should be zero: %2.8f \n',max(abs(StationaryDist0(:)-StationaryDist1(:))))


%%
output=struct(); % Not currently used for anything. Maybe will do so later.

end