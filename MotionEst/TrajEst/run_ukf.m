function outstruct = run_ukf(camstruct, options)
groups = options.est.groups;
link = options.link;
%determine 
%npts = length(pts);
cams = options.cams;
ncam = length(cams);
tstart = options.est.tstart;
tstop  = options.est.tstop;
nsteps = tstop-tstart+options.interp;
%time = linspace(tstart,tstop,nsteps+1);
%Pi0 = [1, 0, 0, 0; 0, 1, 0, 0];

meas = create_meas_matrix(camstruct, options);
%StereoTriangulation_meas_mat(meas,camstruct,options);
%% Filter - EKF

%Covariance Matrix for each point
link_list = get_group_links(options.link, groups);
Sigma_k = calc_Rt_joint(link_list, link);

%State matrix with all points
x_km1 = options.est.state_init;
% x_km2 = x_km1;

% Motion model: Just kinematic point-mass model.  xn = x + u*delta_t;
state_update_model = @(x_km1, x_km2, links, Rt_handle) mm_ConstVel(x_km1, x_km2, options, links, Rt_handle);
%state_jac          = @(x, u, t) dg_ZeroDyn(x,u,t);
Rt_handle          = @(ll)        calc_Rt_joint(ll,options.link);

% Measurement model: Each measurment is a col-vector of length 3*n where n
%   is the number of fixed markers on the bat body.  Assume exact knowledge
%   of the marker positions in body-frame (comes from param.mkr) and exact
%   orientation of the body in 3D space (comes from the upper-left of the H
%   matrix at every timestep)

msmt_model = @(x, t) CamNet_JS(x, camstruct, options);
%msmt_jac   = @(x, t) dh_CamNet(x, camstruct);

%[X_ekf, Sig_X_ekf, mat_mags_ekf] = run_extended_kf(x_km1, Sigma_k, zeros(3*npts,nsteps), meas, state_update_model, state_jac, msmt_model, msmt_jac, Rt_handle);
[X_ukf, Sig_X_ukf]               = run_unscented_kf_recursive2(x_km1, Sigma_k, zeros(options.nstate,nsteps), meas, state_update_model, msmt_model, Rt_handle, options);

%outstruct.ekf.X         = X_ekf;
outstruct.ukf.X         = X_ukf;
%outstruct.ekf.Sig_X     = Sig_X_ekf;
outstruct.ukf.Sig_X     = Sig_X_ukf;

%outstruct.ekf.mat_mags = mat_mags_ekf;
