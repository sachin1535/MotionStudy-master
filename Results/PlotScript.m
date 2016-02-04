function PlotScript(camstruct, inertstruct, options)
tstart      = options.plot.tstart;
tstop       = options.plot.tstop;
time        = linspace(tstart,tstop,tstop-tstart+1);
pts         = options.plot.pts;
%pt_names    = options.pt_names;
npts        = length(pts);
cams        = options.cams;
ncam        = length(cams);
nsteps      = tstop-tstart+1;
linestyle1  = options.plot.linestyle1;
linestyle2  = options.plot.linestyle2;
plot_options = options.plot.fig_txt_props;

%% Plot norm of Sigma_X to Evaluate Convergence
% norm_Sig_X_ekf = zeros(nsteps,1);
% norm_Sig_X_ukf = zeros(nsteps,1);
% 
% for kk = 1:nsteps
%     for pp = 1:npts
%         norm_Sig_X_ekf(kk,pp) = norm(inertstruct.filt.ekf.Sig_X(3*(pp-1)+1:3*pp,3*(pp-1)+1:3*pp,kk));
%         norm_Sig_X_ukf(kk,pp) = norm(inertstruct.filt.ukf.Sig_X(3*(pp-1)+1:3*pp,3*(pp-1)+1:3*pp,kk));
%     end
% end
% 
% figure
% hold on
% for pp = 1:npts
%     plot(time', norm_Sig_X_ekf(:,pp), '-o', time', norm_Sig_X_ukf(:,pp), '-o');
% end
% legend('norm SigX EKF','norm SigBar EKF','norm SigX UKF');

%% Plot Norms of other matricies to determine convergence
% mat_mags.Hj_mag = Hj_mag;
% mat_mags.Qt_mag = Qt_mag;
% mat_mags.Gj_mag = Gj_mag;
% mat_mags.Kk_mag = Kk_mag;
% mat_mags.Knum_mag = Knum_mag;
% mat_mags.Kden_mag = Kden_mag;
% mat_mags.Sig_bar_mag = Sig_bar_mag;

% figure
% plot(time', mat_mags.Hj_mag, '-o', time', mat_mags.Qt_mag, '-o', time', mat_mags.Gj_mag, '-o', time', mat_mags.Kk_mag, '-o');
% legend('||H_j||','||Q_t||', '||G_j||', '||K_k||');
% 
% figure
% plot(time', mat_mags.Kden_mag, '-o')
% legend('||den(K_k)||');
% 
% figure
% plot(time', mat_mags.Knum_mag, '-o')
% legend('||num(K_k)||');
% 
% figure 
% plot(time', mat_mags.KkHj_mag, '-o')
% legend('||K_kH_j||');

%% Plot Error between KF results and Stereo Triangulation.
% error_ekf = zeros(nsteps,npts);
% error_ukf = zeros(nsteps,npts);
% for kk = 1:nsteps
%     for pp = 1:npts
%         error_ekf(kk,pp) = norm(inertstruct.filt.ekf.X(3*(pp-1)+1:3*pp,kk) - inertstruct.ster(1).pts(:,kk,pp)); 
%         error_ukf(kk,pp) = norm(inertstruct.filt.ukf.X(3*(pp-1)+1:3*pp,kk) - inertstruct.ster(1).pts(:,kk,pp));
%     end
% end
% 
% %EKF error
% figure
% hold on
% 
% for pp = 1:npts
%     plot(time',error_ekf(:,pp), linestyle1{pp},time',error_ukf(:,pp),linestyle1{pp}, 'LineWidth', 2)   
% end
% xlabel('Timestep', 'FontSize', 14)
% ylabel('Magnitude of error (mm)', 'FontSize', 14)
% title ('Error: EKF vs Stereo Triagulation', 'FontSize', 14)
% V = axis;
% axis([time(1),time(nsteps),V(3),V(4)])
% set(gca, 'FontSize', 12);
% 
% %UKF Error
% figure
% hold on
% for pp = 1:npts
%     plot(time',error_ukf(:,pp),linestyle1{pp}, 'LineWidth', 2)%,time',error_ukf(:,pp),linestyle_ster{pp})   
% end
% xlabel('Timestep', 'FontSize', 14)
% ylabel('Magnitude of error (mm)', 'FontSize', 14)
% title ('Error: UKF vs Stereo Triagulation', 'FontSize', 14)
% V = axis;
% axis([time(1),time(nsteps),V(3),V(4)])
% set(gca, 'FontSize', 12);

%plot points to visualize
figure
hold on
skel = [];
skel_step = 5;
for pp = 1:npts
    plot3(inertstruct.filt.ekf.X(3*(pp-1)+1,tstart:tstop)', inertstruct.filt.ekf.X(3*(pp-1)+2,tstart:tstop)', inertstruct.filt.ekf.X(3*pp,tstart:tstop)',linestyle1{pp})
    plot3(inertstruct.filt.ukf.X(3*(pp-1)+1,tstart:tstop)', inertstruct.filt.ukf.X(3*(pp-1)+2,tstart:tstop)', inertstruct.filt.ukf.X(3*pp,tstart:tstop)',linestyle2{pp})
%    text(inertstruct.filt.ukf.X(3*(pp-1)+1,5)', inertstruct.filt.ukf.X(3*(pp-1)+2,5)', inertstruct.filt.ukf.X(3*pp,5)-100, pt_names{pp}, plot_options{:}, 'VerticalAlignment', 'bottom')
%     if pts(pp)>12
%         skel = [skel;inertstruct.filt.ekf.X(3*(pp-1)+1:3*pp,skel_step)';inertstruct.filt.ekf.X((3*(find(pts==12)-1)+1):3*find(pts==12),skel_step)'];
%     else
%         skel = [skel;inertstruct.filt.ekf.X(3*(pp-1)+1:3*pp,skel_step)'];
%     end
end

%plot3(skel(:,1),skel(:,2),skel(:,3),'b','LineWidth',2);
handle = legend('EKF', 'UKF');%, 'Stereo 12', 'Stereo 13', 'Stereo 23')
set(handle, 'FontSize', 16);

H = zeros(4,4,ncam);
for cc = 1:ncam
    H(:,:,cc) = camstruct(cams(cc)).H;
end
% CFPlot(H, 50)
axis equal
set(gca, 'FontSize', 16, 'CameraPosition', [0, 0, 0])
xlabel('x (mm)', 'FontSize', 16)
ylabel('y (mm)', 'FontSize', 16)
zlabel('z (mm)', 'FontSize', 16)
title ('Marker Trajectory in 3D Coordinates', 'FontSize', 16)

