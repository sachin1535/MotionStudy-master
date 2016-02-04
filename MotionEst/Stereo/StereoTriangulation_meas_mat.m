function StereoTriangulation_meas_mat(meas, camstruct, options)
ncam      = length(camstruct);
cams      = 1:ncam;
npts      = size(meas,1)/(2*ncam);
pts       = 1:npts;
timesteps = 1:size(meas,2);
nsteps    = length(timesteps);
linestyle1 = options.plot.linespec1;
plot_start = 1;
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
%stereostruct = struct([]);
for pair = 1:npair
    for pp = 1:npts
        for kk = 1:nsteps
            cam1 = pair_list(pair,1);
            cam2 = pair_list(pair,2);
            pt_c1 = meas(2*(cam1-1)*npts+2*(pp-1)+1:2*(cam1-1)*npts+2*pp,timesteps(kk));
            pt_c2 = meas(2*(cam2-1)*npts+2*(pp-1)+1:2*(cam2-1)*npts+2*pp,timesteps(kk));
            stereostruct(pair).pts(:,kk,pts(pp)) = stertridet2(pt_c1, pt_c2, camstruct(cam1),camstruct(cam2));
            stereostruct(pair).cams = pair_list(pair,:);
        end
    end
end

%% Plot Stereo Triangulations Using Different Pairs of Cameras
figure
hold on
for pair = 1:npair
    for pp = 1:npts
        plot3(stereostruct(pair).pts(1,plot_start:nsteps,pts(pp))', stereostruct(pair).pts(2,plot_start:nsteps,pts(pp))', stereostruct(pair).pts(3,plot_start:nsteps,pts(pp))',linestyle1{pp})
    end
end
%legend('PT 2 Cams 1 and 2', 'PT 2 Cams 1 and 3', 'PT 2 Cams 2 and 3')

%CFPlot(H, 50)
axis equal
set(gca, 'FontSize', 16, 'CameraPosition', [0, 0, 0])
xlabel('x (mm)', 'FontSize', 16)
ylabel('y (mm)', 'FontSize', 16)
zlabel('z (mm)', 'FontSize', 16)
title ('Stereo Triangulation', 'FontSize', 18)