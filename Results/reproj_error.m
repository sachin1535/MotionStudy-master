function reproj_error(camstruct, inertstruct, options)
%proj error flight 3 progrram flow 
%% Import Plot Settings
cams = options.est.cams;
camstruct = camstruct(cams);
link_names = options.link_names;
ncam = length(cams);
pts  = options.plot.pts;
npts = length(pts);
tstart = options.tstart;
tstop = options.tstop;
nsteps = tstop-tstart+options.interp;
filepath = options.path;
colors = options.plot.colors;
save_fig    = options.plot.savefig;
save_reproj = options.plot.saveim_reproj;
save_reproje = options.plot.saveim_reproje;
fig_txt_props = options.plot.fig_txt_props;

plot_stop  = options.plot.tstop;
plot_start = options.plot.tstart;
fig_num = options.plot.reprojframe;

%% Transform inertial points into image space
cnt = 0;
for pp = pts
    cnt = cnt+1;
    X_ukf(3*(cnt-1)+1:3*cnt,:)  = inertstruct.ukf.Features(3*(pp-1)+1:3*pp,3:end);
end

y_k_ukf         = zeros(2*ncam*npts,nsteps);
lambdas_ukf     = zeros(npts,nsteps,ncam);
meas_dist       = zeros(2*ncam*npts,nsteps);

%compile measurements into a matrix
meas = zeros(ncam*npts*2,nsteps);
for c = 1:ncam
    %load([options.path,filesep,'..',filesep,'Calibration_run',filesep,'Intrinsic',filesep,'CalTech',filesep,'Cam',num2str(cams(c)),filesep,'int_cam',num2str(cams(c)),'.mat'],'kc')
    %camstruct(c).dist_c = kc;
    t = options.tstart-camstruct(c).start_frame+1+floor(camstruct(c).sync_del*120):options.tstop-camstruct(c).start_frame+1+floor(camstruct(c).sync_del*120);
    for pp = 1:npts
        meas((c-1)*2*npts+2*(pp-1)+1:(c-1)*2*npts+2*pp,:) = camstruct(c).pts_sync(:,t,options.plot.pts_orig(pp));
    end
end

%Re-inject distortion
for kk = 1:nsteps
    meas_dist(:,kk)                           = CameraDistortion(meas(:,kk),camstruct);
    [y_k_ukf(:,kk), lambdas_ukf(:,kk,:)]      = CamNetDistortion(X_ukf(:,kk),camstruct);
end

%reformat point array
points_ukf = zeros(npts,2,nsteps,ncam);
points_meas = zeros(npts,2,nsteps,ncam);
for c = 1:ncam
    for pp = 1:npts
        for kk = 1:nsteps
            points_meas(kk,:,pp,c)  = meas_dist(2*(c-1)*npts+2*(pp-1)+1:2*(c-1)*npts+2*(pp-1)+2,kk)';
            points_ukf (kk,:,pp,c)  = y_k_ukf(2*(c-1)*npts+2*(pp-1)+1:2*(c-1)*npts+2*(pp-1)+2,kk)';
        end
    end
end

%% Determine Reprojection Error (Reprojected - Measured = Error)
reproj_error_ukf = -points_ukf+points_meas;

%% Determine Mean and STD Dev of Reprojection Errors 
%Determine Reprojection Errors in mm
reproj_error_dist_ukf = zeros(nsteps,2,npts,ncam);
for c = 1:ncam
    for pp = 1:npts
        for kk = 1:nsteps
            reproj_error_dist_ukf(kk,:,pp,c)    = reproj_error_ukf(kk,:,pp,c);
        end
    end
end

%% Determine STD Dev of reprojection errors
%cut out occluded points
AllPts_ukf = [];

for c = 1:ncam
    for pp = 1:npts
        %Grab page of points
        col1_ukf = reproj_error_dist_ukf(plot_start:plot_stop,1,pp,c);
        col2_ukf = reproj_error_dist_ukf(plot_start:plot_stop,2,pp,c);

        %Find Occlusions
        occlusions1_ukf = isnan(col1_ukf);
        occlusions2_ukf = isnan(col2_ukf);

        %Cut occlusions
        col1_ukf(occlusions1_ukf) = [];
        col2_ukf(occlusions2_ukf) = [];
        
        %store in3 matrix
        AllPts_ukf = [AllPts_ukf;col1_ukf,col2_ukf];
        
    end
    reproj_mean_ukf(:,:,c) = mean(AllPts_ukf)';
    reproj_std_ukf(:,:,c)  = std(AllPts_ukf)';
    
    AllPts_ukf = [];
end

%Compute Reprojection STD Ellipses
circ_pt_num = 300;
elipse_ukf = zeros(circ_pt_num,2,ncam);
ndev = 1.97;

cam_hands_ukf = zeros(ncam,1);
for c = 1:ncam
    for ii = 1:circ_pt_num
        elipse_ukf(ii,:,c) = (reproj_mean_ukf(:,:,c) + ndev*reproj_std_ukf(:,:,c).*[cos(2*pi*ii/circ_pt_num);sin(2*pi*ii/circ_pt_num)])';
    end
end

%% Create Plots
x = 1;
y = 2;
reprojerror_fig_ids = zeros(ncam,1);
reproj_fig_ids = zeros(ncam,1);
for c = 1:ncam
    %Project Identified Trajectories onto Images
    %figure
    %If the video has been delaced
    path_name = strcat(options.path,filesep,'Cam',num2str(cams(c)));
    if exist(path_name,'file')
          
        reproj_fig_ids(c) = figure;
        imshow([path_name,filesep,num2str(fig_num+round(camstruct(c).sync_del*119.88)),'.png'])
        hold on
        %Plot measured, EKF, and UKF points
        step = nsteps-5;
        cnt = 0;
        for pt = 1:npts %[1,2,3,6,7,8,11,12,14,15,17];
            cnt = cnt+1;
            p1 = plot(points_meas(plot_start:plot_stop,x,pt,c),points_meas(plot_start:plot_stop,y,pt,c),'+','Color',options.plot.colors2(cnt,:), 'LineWidth', 1);
            p3 = plot(points_ukf(plot_start:plot_stop,x,pt,c),points_ukf(plot_start:plot_stop,y,pt,c),'-','Color',options.plot.colors2(cnt,:), 'LineWidth', 1);
            if cnt == 1
                legend_handles = [p1 p3];
            end
        end
        %[hleg1, hobj1] = legend(legend_handles, 'Image Features', 'UKF Reprojection','Location','NorthWest');
        %textobj = findobj(hobj1, 'type', 'text');
        %lineobj = findobj(hobj1, 'type', 'line');
        %set(lineobj, 'LineWidth', 2);
        %set(textobj, fig_txt_props{:});
        %plot settings **These settings are created for 7_16_14 data run 2
        %functionality will be added so these setting can be changed.
        
        if cams(c) == 301
            limits = [550 880 400 720];
        elseif cams(c) ==302
            limits = [500 880 300 525];
        elseif cams(c) ==303
            limits = [413 890 316 650];
        elseif cams(c) ==310
            limits = [400 900 150 450];    
        elseif cams(c) ==312
            limits = [100 1000 280 720];
        elseif cams(c) ==318
            limits = [450 1050 100 400];
        elseif cams(c) ==320
            limits = [40 900 340 670];
        elseif cams(c) ==325
            limits = [519 880 175 420];
        elseif cams(c) ==333
            limits = [325 570 275 480];    
        else
            limits = [0 1280 0 720];
        end
        
        title(['Cam ',num2str(c)], fig_txt_props{:})
        %set(handle, 'Color', 'w', 'FontSize', 20)
        axis(limits)
        frame_struct = getframe;
        
        %Do the figures need to be saved?
        if save_reproj
            imwrite(frame_struct.cdata,strcat(options.plot.savepath,filesep,'BatFlight_EKFUKF_Reproj_c',num2str(c),'.png'),'png')
        end
        if save_fig
            saveas(reproj_fig_ids(c),strcat(options.plot.savepath,filesep,'BatFlight_EKFUKF_Reproj_c',num2str(c),'.fig'),'fig')
        end

    end
    reprojerror_fig_ids(c) = figure;
    hold on
    cnt = 0;
    for pp = 1:npts%[1,2,3,6,7,8,11,12,14,15,17]
        cnt = cnt+1;
        p2 = plot(reproj_error_dist_ukf(plot_start:plot_stop,1,pp,c),reproj_error_dist_ukf(plot_start:plot_stop,2,pp,c),'o','Color',options.plot.colors2(cnt,:));
        if cnt == 1
            legend_handles = p2;
        end
    end
    
    p4 = plot(elipse_ukf(:,1,c),elipse_ukf(:,2,c), '--k', 'LineWidth', 1.5);
    axis equal
    h(c) = gca;
    xlims(c,:) = xlim;
    ylims(c,:) = ylim;
    
    set(h(c),fig_txt_props{:})
    textobj = findobj(h(c), 'type', 'text');
    lineobj = findobj(h(c), 'type', 'line');
    set(lineobj, 'LineWidth', 2);
    set(textobj,  fig_txt_props{:});
    
    %handle = legend(legend_handles, 'UKF', '95% Conf', 'Location','NorthEast');
    %set(handle,  fig_txt_props{:});
    xlabel('x (pixels)',fig_txt_props{:}); ylabel('y (pixels)',  fig_txt_props{:}); 
    title(sprintf('Cam %d',c));
    %set(textobj, 'Interpreter', 'latex', 'fontsize', 20, 'FontWeight', 'Bold');
end

x_min = min(min(xlims([1,3,8,9],:)));
x_max = max(max(xlims([1,3,8,9],:)));
y_min = min(min(ylims([1,3,8,9],:)));
y_max = max(max(ylims([1,3,8,9],:)));

axis(h([1,3,8,9]), [x_min, x_max, y_min, y_max]);

for im = 1:length(h)
    set(reprojerror_fig_ids(im),'color','w');
    frame_struct = getframe(reprojerror_fig_ids(im));
    if save_reproje
        %imwrite(frame_struct.cdata,strcat(options.plot.savepath,'BatFlight_EKFUKF_ReprojError_c',num2str(im),'.png'),'PNG')
        print('-r600', '-dpdf', strcat(options.plot.savepath,filesep,'BatFlight_EKFUKF_ReprojError_c',num2str(im),'.pdf'))
    end
    if save_fig
        saveas(reprojerror_fig_ids(c),strcat(options.plot.savepath,filesep,'BatFlight_EKFUKF_ReprojError_c',num2str(im),'.fig'),'fig')
    end
end


    