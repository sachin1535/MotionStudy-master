%clear all
%close all
%% --------------------------Set Options----------------------------------
%Import Options
options.groups          = [1,2,3];
options.link_names      = {'Body','Humerus','Radius','Metacarpal 3', 'Metacarpal 4','Metacarpal 5'};
options.tstart          = 390;                  %Note: due to sync delay the first 
options.tstop           = 425;                  %Useable timestep will be tstart+1 
options.interp          = 1;                    %1- data Was NOT interpolated, 0- otherwise;
options.cams            = [301,302,303,310,312,318,320,325,333];
options.plotflag        = 0;
options.path            = 'C:\ShandongData2015\Batflight_07242015\Test004';
options.default_dir     = pwd;

%Bundle Adjustment Options
% options.ba                  = 1;                          %3 load previous and run BA, 2 run BA from scratch, 1 load previous, 0 no BA
% options.ba_tsteps           = [3, 8, 17, 21];             %Row vec of tsteps to use in BA
% options.ba_pts              = [2, 3, 10, 16, 20, 23];     %Row vec of pts to use in BA

%Stereo Options
% options.stereo.links        = options.links;
% options.stereo.cams         = options.cams;
% options.stereo.tstart       = 1;
% options.stereo.tstop        = options.tstop-options.tstart;

%Trajectory Estimation Options
options.est.groups          = options.groups;
options.est.tstart          = 1;
options.est.tstop           = options.tstop - options.tstart+1;
options.est.state_init      = [-1.097,-0.2412,-0.2512,-pi/2,150*pi/180,pi,...%]';%,...
                                0,-pi/2,pi/2,...
                                30/180*pi...
                                pi/4,-20*pi/180,...
                                pi/2,-pi/4,...
                                100/180*pi,-20*pi/180]';

%Plot Options
options.plot.pts           = [1:16];
options.plot.reprojframe   = 400;
options.plot.tstart        = 5;
options.plot.tstop         = (options.tstop - options.tstart)-(options.plot.tstart-1);
options.plot.linespec1        = {'.-r','.-b','.-g', '.-m','.-k','.-c','.--r','.--b','.--g','+-r','+-b','+-g', '+-m','+-k','+-c','+--r','+--b','+--g'};
options.plot.linespec2        = {'+-r','+-b','+-g', '+-m','+-k','+-c','+--r','+--b','+--g','o-r','o-b','o-g', 'o-m','o-k','o-c','o--r','o--b','o--g'};
options.plot.linespec3        = {'o-r','o-b','o-g', 'o-m','o-k','o-c','o--r','o--b','o--g','.-r','.-b','.-g', '.-m','.-k','.-c','.--r','.--b','.--g'};
options.plot.colors         = {'r', 'g', 'b', 'c', 'm', 'k'};
options.plot.savepath       = 'D:\Users\Matt\Documents\VT\Research\Motion_estimation_dev\Papers_and_Presentations\Bender2016SciTechPres\';
options.plot.savefig        = 1;
options.plot.saveim_reproj  = 1;
options.plot.saveim_reproje = 1;
options.plot.fig_txt_props  = {'FontName', 'Times New Roman', 'FontSize', 18, 'FontWeight', 'Bold'};

%% Define the Skeleton
SkeletonDefn_BatFlight_20150724_test004
links        = get_group_links(synthConfig.link,options.groups);
options.link = synthConfig.link(links);
options      = create_state_vec(options);
options      = create_meas_vec(options);

%% Load The Camera Measurements
load([options.path,filesep,'CamStruct.mat'])
Cam = Cam(options.cams);
ncam = length(Cam);

for cc = 1:ncam
    dt = Cam(cc).start_frame-1+floor(119.88*Cam(cc).sync_del);
    Cam(cc).pts = [];
    Cam(cc).pts(:,:,1:3) = Cam(cc).pts_sync(:,options.tstart-Cam(cc).start_frame+1+floor(Cam(cc).sync_del*119.88):options.tstop-Cam(cc).start_frame+1+floor(Cam(cc).sync_del*119.88),[1,4,5]);
    Cam(cc).pts(:,:,4:5) = Cam(cc).pts_sync(:,options.tstart-Cam(cc).start_frame+1+floor(Cam(cc).sync_del*119.88):options.tstop-Cam(cc).start_frame+1+floor(Cam(cc).sync_del*119.88),[6:-1:5]);
    Cam(cc).pts(:,:,6:8) = Cam(cc).pts_sync(:,options.tstart-Cam(cc).start_frame+1+floor(Cam(cc).sync_del*119.88):options.tstop-Cam(cc).start_frame+1+floor(Cam(cc).sync_del*119.88),[8:-1:6]);
    
    Cam(cc).pts(:,:,9:11) = Cam(cc).pts_sync(:,options.tstart-Cam(cc).start_frame+1+floor(Cam(cc).sync_del*119.88):options.tstop-Cam(cc).start_frame+1+floor(Cam(cc).sync_del*119.88),[10:-1:8]);
    Cam(cc).pts(:,:,12:14) = Cam(cc).pts_sync(:,options.tstart-Cam(cc).start_frame+1+floor(Cam(cc).sync_del*119.88):options.tstop-Cam(cc).start_frame+1+floor(Cam(cc).sync_del*119.88),[14,13,8]);
    Cam(cc).pts(:,:,15:16) = Cam(cc).pts_sync(:,options.tstart-Cam(cc).start_frame+1+floor(Cam(cc).sync_del*119.88):options.tstop-Cam(cc).start_frame+1+floor(Cam(cc).sync_del*119.88),[17,8]);
    Cam(cc).pt_assoc = [1,1,1,2,2,3,3,3,4,4,4,5,5,5,6,6;1,2,3,1,2,1,2,3,1,2,3,1,2,3,1,2];
end



