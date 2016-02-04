%clear all
close all
%% --------------------------Set Options----------------------------------
%Import Options
options.groups          = [1,2,3];
options.link_names      = {'body', 'Right Humerus', 'Right Radius', 'Metacarpal 1', 'Metacarpal 2','Metacarpal 3' };
options.tstart          = 1;                  %Note: due to sync delay the first 
options.tstop           = 100;                  %Useable timestep will be tstart+1 
options.interp          = 1;                    %1- data Was NOT interpolated, 0- otherwise;
options.cams            = [1:40];
options.plotflag        = 0;
options.path            = 'C:\ShandongData2015\Batflight_07242015\Synthetic';
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

%Plot Options
% options.plot.links           = options.links;
% options.plot.reprojframe    = 220;
% options.plot.tstart         = 5;
% options.plot.tstop          = (options.tstop - options.tstart)-(options.plot.tstart-1);
options.plot.linespec1        = {'.-r','.-b','.-g', '.-m','.-k','.-c','.--r','.--b','.--g','+-r','+-b','+-g', '+-m','+-k','+-c','+--r','+--b','+--g'};
options.plot.linespec2        = {'+-r','+-b','+-g', '+-m','+-k','+-c','+--r','+--b','+--g','o-r','o-b','o-g', 'o-m','o-k','o-c','o--r','o--b','o--g'};
options.plot.linespec3        = {'o-r','o-b','o-g', 'o-m','o-k','o-c','o--r','o--b','o--g','.-r','.-b','.-g', '.-m','.-k','.-c','.--r','.--b','.--g'};
% options.plot.colors         = {'r', 'g', 'b', 'c', 'm', 'k'};
% options.plot.savepath       = 'D:\Users\Matt\Documents\VT\Research\Motion_estimation_dev\Papers_and_Presentations\Bender2015scitech\';
% options.plot.savefig        = 0;
% options.plot.saveim_reproj  = 0;
% options.plot.saveim_reproje = 0;
% options.plot.fig_txt_props  = {'FontName', 'Times New Roman', 'FontSize', 18, 'FontWeight', 'Bold'};

%% Define the Skeleton
SkeletonDefn2
links = get_group_links(synthConfig.link,options.groups);
options.link = synthConfig.link(links);
%determine State Indicies
options = create_state_vec(options);
options = create_meas_vec(options);

load([options.path,filesep,'CamStruct.mat'])
Cam = Cam(options.cams);
ncam = length(Cam);
for cc = 1:ncam
    Cam(cc).pts = Cam(cc).pts(:,options.tstart:options.tstop,1:length([options.link(links).BFvecs]));
    Cam(cc).pt_assoc = Cam(cc).pt_assoc(:,ismember(Cam(cc).pt_assoc(1,:),links));
end
