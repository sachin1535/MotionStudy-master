function stereostruct = StereoTriangulation_svob(camstruct, options)
ncam      = length(options.stereo.cams);
cams      = options.stereo.cams;
npts      = length(options.stereo.pts);
pts       = options.stereo.pts;
timesteps = options.stereo.tstart:options.stereo.dt:options.stereo.tstop;
nsteps    = length(timesteps);
linestyle1 = options.plot.linestyle1;
plot_start = options.stereo.tstart;
fs = 120;
%% Perform Stereo Triangulation for Comparision
npair = 0;
pair_list = [];
for ii = 1:ncam-1
    npair = npair + ii;
    for jj = ii+1:ncam
        pair_list = [pair_list; cams(ii),cams(jj)];
    end
end

%Determine stereo reconstructions of each pair of cams with new extrinsics
stereostruct = struct([]);
for pair = [1:npair]
    for pp = 1:npts
        for kk = 1:nsteps
            %if ~isempty(intersect(pp,camstruct(pair_list(pair,1)).idin(1,pp))) && ~isempty(intersect(pp,camstruct(pair_list(pair,2)).idin(1,pp)))
                camstruct1 = camstruct(pair_list(pair,1));
                camstruct2 = camstruct(pair_list(pair,2));
                stereostruct(pair).pts(:,kk,pts(pp)) = stertridet2(camstruct1.pts_sync(:,timesteps(kk)-camstruct1.start_frame+1+floor(camstruct1.sync_del*fs),pts(pp)), ...
                                                   camstruct2.pts_sync(:,timesteps(kk)-camstruct2.start_frame+1+floor(camstruct2.sync_del*fs),pts(pp)), ...
                                                   camstruct1,camstruct2);
            
            stereostruct(pair).cams = pair_list(pair,:);
        end
    end
end

%% Plot Stereo Triangulations Using Different Pairs of Cameras
figure
hold on
for pair = 1:npair
    if ~isempty(stereostruct(pair).pts)
        for pp = 1:npts
            %plot3(stereostruct(pair).pts(1,options.ba_tsteps,pts(pp))', stereostruct(pair).pts(2,options.ba_tsteps,pts(pp))', stereostruct(pair).pts(3,options.ba_tsteps,pts(pp))','.r')
            plot3(stereostruct(pair).pts(1,:,pts(pp))', stereostruct(pair).pts(2,:,pts(pp))', stereostruct(pair).pts(3,:,pts(pp))',linestyle1{pp})
        end
    end
end
%legend('PT 2 Cams 1 and 2', 'PT 2 Cams 1 and 3', 'PT 2 Cams 2 and 3')
H = reshape([camstruct(cams).H],4,4,[]);
CFPlot(H, 0.1)
axis equal
set(gca, 'FontSize', 16, 'CameraPosition', [0, 0, 0])
xlabel('x (mm)', 'FontSize', 16)
ylabel('y (mm)', 'FontSize', 16)
zlabel('z (mm)', 'FontSize', 16)
title ('Stereo Triangulation', 'FontSize', 18)

%----------------------------STERTRIDENT2-----------------------------
function x = stertridet2(phi1, phi2, camstruct1,camstruct2)

H1 = camstruct1.H;
H2 = camstruct2.H;
K1 = camstruct1.K;
K2 = camstruct2.K;
%Determine Vector from o1 to o2 in world frame.
P1= H1(1:3,4);
P2 = H2(1:3,4);
T = P2-P1;

%determine unit vectors to observed by cameras in world frame
 u1 = H1(1:3,1:3)*(K1\[phi1;1])/norm((K1\[phi1;1]));
 u2 = H2(1:3,1:3)*(K2\[phi2;1])/norm((K2\[phi2;1]));

%Determine Plane parallel to both lines. 
N_hat = cross(u1,u2);

%Project vector between coordinate frames onto normal vector (i.e. the distance between planes)
d_vec = T'*N_hat*N_hat;

T_til = T-d_vec;

gamma = acos(u1'*u2/(norm(u1)*norm(u2)));
beta  = pi-acos(T_til'*u2/(norm(T_til)*norm(u2)));
alpha = acos(T_til'*u1/(norm(T_til)*norm(u1)));
b = norm(T_til)/sin(gamma)*sin(beta)*u1+P1;
a = norm(T_til)/sin(gamma)*sin(alpha)*u2+P2;

%d_vec_new = b-(a+T);
x = b/2+a/2;


