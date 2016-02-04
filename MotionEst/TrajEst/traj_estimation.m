function [eststruct, options] = traj_estimation(camstruct, options)
%% Run Initialization Script for Desired Data Set
Init_Flight_BatFlight_20150724_test004_wrist

%StereoTriangulation_svob(Cam,options);
%% Run UKF 
fprintf('Trajectory Estimation is Running .... \n')
tic
eststruct = run_ukf(camstruct, options);
toc
fprintf('Trajectory Estimation is Complete .... \n')

clear title xlabel ylabel
q = eststruct.ukf.X;
nsteps = size(q,2);
%create the body basis transforms and the coordinate transforms for each
%time step2
links = get_group_links(options.link,options.groups);
dof = 0;
for ll = links
    nDof = options.link(ll).nDof;
    dof = dof+options.link(ll).nDof;
    for kk = 1:nsteps
        q_lk = q(dof-nDof+1:dof,kk);
        if strcmp(options.link(ll).IDkern,'DH')
            kinc(kk).link(ll).H=DHTransforms(q_lk+options.link(ll).thetas,...
                        options.link(ll).alphas,...
                        options.link(ll).disps,...
                        options.link(ll).offsets);
        elseif  strcmp(options.link(ll).IDkern,'YPR')
            theta = q_lk(4:6);
            d = q_lk(1:3);
            kinc(kk).link(ll).H = YPRTransform(theta, d);
        end
    end
end

npts = size([options.link.BFvecs],2);
features = zeros(3*npts,nsteps);
n_pts_tot = 0;
for ll = links
    n_bf_pts = size(options.link(ll).BFvecs,2);
    n_pts_tot = n_pts_tot+n_bf_pts;
    for kk = 1:nsteps
    X = [eye(3,3),zeros(3,1)]*hnode2node(kinc(kk),options,1,ll)*[options.link(ll).BFvecs;ones(1,n_bf_pts)];
    features(3*(n_pts_tot-n_bf_pts)+1:3*n_pts_tot,kk) = X(:);
    end
    assoc(:,(n_pts_tot-n_bf_pts)+1:n_pts_tot) = [ll*ones(1,n_bf_pts);1:n_bf_pts];
end

eststruct.ukf.Features = features;
eststruct.kinc = kinc;