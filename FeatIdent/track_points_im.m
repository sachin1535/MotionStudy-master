function camstruct = track_points_im(camstruct,path)
cam_num = input('Enter the CAMERA number you wish to track a point in: ');
Cam = camstruct(cam_num);

pts = input('Enter the POINT numbers you wish to track: ');

tstart = input('Enter the TIMESTEP to START tracking:');
tstop  = input('Enter the TIMESTEP to STOP tracking:');%Cam.end_frame;
  
%% run the filter
[out, Cam] = run_ukf(Cam,tstart,tstop,pts);

%% plot the results
figure
hold on
for pp = 1:length(pts)
    plot(out.ukf.X(2*pp-1,:),out.ukf.X(2*pp,:))
    plot(Cam.pts(1,:,pts(pp)),Cam.pts(2,:,pts(pp)),'or')
end

camstruct(cam_num) = Cam;


%--------------------------------RUN_UKF-----------------------------------
function [outstruct, camstruct] = run_ukf(camstruct, tstart,tstop,pts)
npts = length(pts);
%cams = options.cams;
%ncam = length(cams);
nsteps = tstop-tstart+1;
options.tstop = tstop;
options.tstart = tstart;
options.path = 'C:\ShandongData2015\Batflight_07242015\Test004';

meas = zeros(npts*2,nsteps);
steps = tstart-camstruct.start_frame+1:tstart-camstruct.start_frame+2;
for pp = 1:npts
    for kk = steps
        meas(2*(pp-1)+1:2*pp,kk-steps(1)+1) = camstruct.pts(:,kk,pts(pp));
    end
end

%seed current covariance 
Sigma_k = eye(2*npts);

% Motion model: Just kinematic point-mass model.  xn = x + u*delta_t;
state_update_model = @(x, t) mm_spline(x,t);
Rt_handle          = @()        calc_Rt_pix();

% Measurement model: Each measurment is a col-vector of length 3*n where n
%   is the number of fixed markers on the bat body.  Assume exact knowledge
%   of the marker positions in body-frame (comes from param.mkr) and exact
%   orientation of the body in 3D space (comes from the upper-left of the H
%   matrix at every timestep)

msmt_model = @(x, kk) prox_model(x, kk, camstruct, options);

[X_ukf, Sig_X_ukf, z]        = run_unscented_kf(meas(:,2), Sigma_k, zeros(3*npts,nsteps), meas, state_update_model, msmt_model, Rt_handle);

%outstruct.ekf.X         = X_ekf;
outstruct.ukf.X         = X_ukf;
%outstruct.ekf.Sig_X     = Sig_X_ekf;
outstruct.ukf.Sig_X     = Sig_X_ukf;
%camstruct.pts = zeros(2,size(z,2));
for pp = 1:length(pts)
    camstruct.pts(:,tstart-camstruct.start_frame+1:tstop-camstruct.start_frame+1,pts(pp)) = z(2*pp-1:2*pp,:);
end
%outstruct.ekf.mat_mags = mat_mags_ekf;


%----------------------------RUN_UNSCENTED_KF------------------------------
function [X, Sig_X, z, varargout] = run_unscented_kf(...
    mu_0, Sig_0, u, z, g_handle, h_handle, Rt_handle, varargin)
% RUN_UNSCENTED_KF
%  [X, Sig_X] = RUN_UNSCENTED_KF(mu_0, Sig_0, u, z, g_handle, h_handle)
%    runs unscented kalman filter on input datasets u and z (described
%    below) with input motion and sensor models (also described below)
%    
%    INPUTS:
%    mu_0 = Initial state estimate vector (size nx1)
%    Sig_0 = Initial covariance estimate (square pos-def matrix size nxn)
%    u = Input trajectory matrix (size mxT, T = num-of-timesteps)
%    z = Msmt matrix (size pxT, T = num-of-timesteps)
%    g_handle = func handle to Motion model, see below for form
%    h_handle = func handle to Sensor model, see below for form
%
%    OUTPUTS:
%    X = State estimate matrix (size nxT) for all time T
%    Sig_X = State covariance (size nxnxT) for all time T
%
%    FUNCTION HANDLES:
%    [xk, Rk] = g_handle(x,u,t) returns the next state prediction (xk) from
%      some prediction model on the current state x, the input u, and
%      (optionally) the timestep t.  The dimensions of xk must be equal to
%      the dimensions of x (also equivalent to dim mu_0 above).  Also
%      returns the associated covariance Rk.
%    [zk, Qk] = h_handle(x,t) returns the measurement prediction (zk) from
%      some measurement model on the current state x and (optionally) the
%      timestep t.  Also returns the associated covariance Qk.
% 
%  [..., outstruct] = RUN_UNSCENTED_KF(..., param) allows the input of
%    an optional param struct to control the algorithm and the output of
%    optional params the user may want to investigate performance
%
%    param.lam = parameter controling spread of sigma pts.  lam=1 by
%     default, (lam+n) = 3 good for gaussians (according to julier1997)
%
%    outstruct.mu_bar = [nxT] matrix of motion-model state estimates
%    outstruct.Sig_bar = [nxnxT] matrix of motion-model state covariances 

if nargin < 8
    param = struct();
else
    param = varargin{1};
end

path = 'C:\ShandongData2015\Batflight_07242015\Test004';
plot_flag = 0;

% Note: code follows variable naming in thrun2005probab_robot textbook

% Specify default parameters
default.lam = 1; % Controls spread of sigma points

default.beta = 2; % Advanced param, optimal = 2 for gaussians
default.alpha = .25; %1; % Advanced param, value guessed by hgm

% Use passed-in params, or defaults if blank
param = populate_struct_with_defaults(param, default);

n = length(mu_0);
m = 2*n+1;

% build w col-vectors (weights of sigma points)
wm(1) = param.lam/(n + param.lam);
wc(1) = wm(1) + (1-param.alpha^2 + param.beta);
wm(2:2*n+1) = 1/(2*(n+param.lam))*ones(2*n,1);
wc(2:2*n+1) = wm(2:2*n+1);

% Initialize outputs
X = zeros(n,size(z,2));
X(:,1:2) = z(:,1:2);
Sig_X = zeros(n,n,size(z,2));

mu = mu_0;
Sig = Sig_0;
% norm_Q_log = zeros(size(u,2),1);

for ii = 3:size(u,2) % for all timesteps
    % Line 2: Unscented transform on mu, Sig
    if ii ==1
        niter = 30;
    else
        niter = 1;
    end
    
    for zz = 1:niter
        if plot_flag
            figure(1)
            clf;
            imshow(imread([path,filesep,'Cam309',filesep,num2str(371+ii),'.png']))
            hold on
        end
        Chi_prev = unscented_transform(mu, param.lam, Sig);
        
        % Line 3: Form Chi_star via motion model on sig points
        Chi_star = zeros(size(Chi_prev));
        for sigpt = 1:size(Chi_star, 2)
            [Chi_star(:,sigpt), ~] = g_handle([X(:,1:ii-2),Chi_prev(:,sigpt)], [1:ii-1]);
        end

        % Line 4&5: Form mu_bar and Sig_bar from weighted sum of Chi_star
        mu_bar = Chi_star * wm';
        if plot_flag
            plot(X(1:2:end,ii-1),X(2:2:end,ii-1),'xy')
            plot(mu_bar(1:2:end),mu_bar(2:2:end),'+r')
        end
        [~, Rt] = g_handle([X(:,1:ii-2),mu_bar], [1:ii-1]);
        del_Sig = (Chi_star - repmat(mu_bar,1,size(Chi_star,2)));
        Sig_bar = zeros(size(Sig));
        for ndx = 1:size(del_Sig,2)
            Sig_bar = Sig_bar + wc(ndx)*del_Sig(:,ndx)*del_Sig(:,ndx)';
        end
        Sig_bar = Sig_bar + Rt;

        % Line 6: Chi_bar from Unscented transform
        Chi_bar = unscented_transform(mu_bar, param.lam, Sig_bar);

        % Line 7: Z_bar from msmt model on sig points
        Z_bar = zeros(size(z,1), size(Chi_bar,2));
        for sigpt = 1:size(Chi_bar, 2)
            Z_bar(:,sigpt) = eye(n)*Chi_bar(:,sigpt);
        end
        
        % Detect & handle occlusions
        z(:,ii) = h_handle(Z_bar(:,1),ii);
        occlusion_ndx = find(isnan(z(:,ii))); % find ndx of occlusions
        z_minus_occlusions = z(:,ii); % create local msmt copy
        z_minus_occlusions(occlusion_ndx) = []; % strip occlusions out
        
        if plot_flag
            plot(z_minus_occlusions(1:2:end),z_minus_occlusions(2:2:end),'oc')
        end
              
        Z_bar(occlusion_ndx,:) = []; % strip out occluded measurments

        % Line 8--10: z_hat (mean msmt) and S (msmt cov) from weighted sum
        % of Z_bar, Sig_xz from "del terms"
        z_hat = Z_bar * wm';
        if plot_flag
            plot(z_hat(1:2:end),z_hat(2:2:end),'sg')
        end
        
        Qt = 10*eye(n);
        del_Z = (Z_bar - repmat(z_hat,1,size(Z_bar,2)));
        del_Sig = Chi_bar - repmat(mu_bar,1,size(Chi_bar,2));
        S = zeros(length(z_hat), length(z_hat));
        Sig_xz = zeros(length(mu_bar), length(z_hat));
        for ndx = 1:size(del_Z,2)
            S = S + wc(ndx)*del_Z(:,ndx)*del_Z(:,ndx)';
            Sig_xz = Sig_xz + wc(ndx)*del_Sig(:,ndx)*del_Z(:,ndx)';
        end

        Qt(occlusion_ndx,:) = []; % strip occlusion rows out
        Qt(:, occlusion_ndx) = []; % strip occlusion cols out
        S = S + Qt;

        % Line 11: K (kalman msmt gain)
        K = ((S')\Sig_xz')'; % just avoid inv(S)

        % Line 12: Update state measurement and cov
        mu = mu_bar + K*(z_minus_occlusions-z_hat); % use stripped msmt
        if plot_flag
            plot(mu(1:2:end),mu(2:2:end),'^m')
            pause
        end
        Sig = Sig_bar - K*S*K';

        % Accumulate measurements into function outputs
        X(:,ii) = mu;
        Sig_X(:,:,ii) = Sig;
        % norm_Q_log(ii) = max(max(Qt));

        outstruct.mu_bar(:,ii) = mu_bar;
        outstruct.Sig_bar(:,:,ii) = Sig_bar;
    end
end
% norm_Q_log
if nargout > 2
    varargout{1} = outstruct;
end
    
%------------------------------UNSCENTED_TRANSFORM-------------------------
function Chi = unscented_transform(mu, lam, Sig)
Chi = zeros(length(mu), 2*length(mu)+1);
Chi(:,1) = mu; % sigpt_0 is Chi_prev(1)
sqrtSig = sqrtm(Sig); % Matrix square-root.  Possibly could do faster
                      % with Cholesky decomp? (for future ref)
for ndx = 1:length(mu)
    sigpt_adj = sqrt(length(mu)+lam)*sqrtSig(:,ndx);
    Chi(:,ndx+1) =            mu + sigpt_adj;
    Chi(:,ndx+1+length(mu)) = mu - sigpt_adj;
end

%------------------------------PROX_MODEL---------------------------------
function [y, Qk] = prox_model(x, kk, camstruct, options)
nstates = length(x);
ff = options.tstart - camstruct.start_frame + kk;
y = zeros(size(x));
for pp = 1:nstates/2
    if isempty(camstruct.features{ff})
        y = NaN*ones(length(x),1);
    else
        locs = camstruct.features{ff}.Location;
        inds = find(camstruct.features{ff}.SignOfLaplacian>0);
        locs(inds,:) = []; 
        dist = (locs'-repmat(x(2*pp-1:2*pp,1),1,size(locs,1)));
        delta = sum(dist.*dist,1);
        [val,I] = min(delta);
        if val>10^2
            y(2*pp-1:2*pp) = [NaN;NaN];
        else
            y(2*pp-1:2*pp) = locs(I,:)';
        end
        
    end
end
% if ~isempty(camstruct.features{ff})
%     figure(5)
%     imshow([options.path,filesep,'Cam309',filesep,num2str(ff+camstruct.start_frame-1),'.png'])
%     hold on
%     title(sprintf('Image %d',ff+camstruct.start_frame-1));
%     plot(camstruct.features{ff}.Location(:,1),camstruct.features{ff}.Location(:,2),'oc')
%     plot(x(1:2:nstates),x(2:2:nstates),'+r')
%     plot(y(1:2:nstates),y(2:2:nstates),'^g')
%     pause
% end

Qk = 10*eye(length(y));

%-----------------------------POPULATE_STRUCT_WITH_DEFAULTS----------------
function s = populate_struct_with_defaults(s, default_struct)
% POPULATE_STRUCT_WITH_DEFAULTS
%   POPULATE_STRUCT_WITH_DEFAULTS(s, default_struct) populates the input
%   struct s.  The code scans the fields of default_struct, and if any
%   field in s is missing or empty, the default_struct field-value is
%   used.
%
%   Example:
%     d.a = 3; d.b = 4; d.c = 5;
%     s.a = 10; s.b = [];
%     p = populate_struct_with_defaults(s,d);
%   Results in:
%     p.a = 10; (value from s.a)
%     p.b = 4; (value from d.b)
%     p.c = 5; (value from d.c, note that field 'c' was created)
%
%   HGM, 2014-10-07
%
fld_lst = fields(default_struct);
for fld_ndx = 1:numel(fld_lst)
    fld_name = fld_lst{fld_ndx};
    if ~isfield(s, fld_name)
        s.(fld_name) = [];
    end
    if isempty(s.(fld_name));
        s.(fld_name) = default_struct.(fld_name);
    end
end

%--------------------------------MM_SPLINE--------------------------------
function [x_k, R_k] = mm_spline(x_kprev, t)
%MM_SPLINE      -Creates a spline prediction of all states within X_KPREV.
%X_PREV should be a [N X K-1] matrix of state values where N is the number 
%of states in the system and K is the time step being predicted.  The
%function creates a spline fit of all data on each dimension and predicts
%each DOF one step forward. T is a [1XK-1 time vector]
% if length(t)>4
%     t = t(end-3:end);
%     x_kprev = x_kprev(:,end-3:end);
% end
% 
[nstates, ~] = size(x_kprev);
% x_k = zeros(nstates,1);
% for nn = 1:nstates/2
%     x_coords = x_kprev(2*(nn-1)+1,:);
%     y_coords = x_kprev(2*nn,:);
%     x_new = 2*x_coords(end)-x_coords(end-1);
%     y_new = 2*y_coords(end)-y_coords(end-1);
%     x_k(2*(nn-1)+1) = (spline(y_coords,x_coords,y_new)+x_new)/2;
%     x_k(2*nn)       = (spline(x_coords,y_coords,x_new)+y_new)/2;
% end

% dx = x_kprev(:,end)-x_kprev(:,end-1);
% avg_dx = [mean(dx(1:2:nstates));mean(dx(2:2:nstates))];
% delta = zeros(nstates,1);
% delta(1:2:nstates,1) = avg_dx(1);
% delta(2:2:nstates,1) = avg_dx(2);
% x_k = x_kprev(:,end) + delta;

x_k = 2*x_kprev(:,end)-x_kprev(:,end-1);

R_k = 50*eye(nstates);

%----------------------------------CALC_RT---------------------------------
function Rt = calc_Rt()
% This function generates the additive uncertainty for the motion model
% x_t = eye()*x_{t-1} for bat-flight points based on some
% back-of-the-envelope assumptions about how far each point could move
% per frame.  The idea is to construct a 95% confidence ellipse for the
% new location of the point based on a rough approximation of it's
% maximum movement in a single frame.

% A number of assumptions are inherent in this model which may be
% revisited or revised, I hope those are evident from the names.  In
% brief, this model assumes rigid-wing flapping and calculates the
% displacement of a point on the tip of the rigid wing, then square-sums
% that with the displacement from the bat-flight through the world.

Rt = eye(10);
