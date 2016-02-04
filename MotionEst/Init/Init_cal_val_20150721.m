%clear all
close all
%% --------------------------Set Options----------------------------------
%Import Options
options.pts             = [1];
options.pt_names        = {'LED'};
options.tstart          = 1;                  %Note: due to sync delay the first 
options.tstop           = 551;                  %Useable timestep will be tstart+1 
options.interp          = 1;                    %1- data Was NOT interpolated, 0- otherwise;
options.cams            = [301:305];
options.plotflag        = 0;
options.path            = 'C:\ShandongData2015\batflight_07242015\Calibration_run\Extrinsic';
options.default_dir     = pwd;

%Bundle Adjustment Options
options.ba                  = 0;                          %3 load previous and run BA, 2 run BA from scratch, 1 load previous, 0 no BA
options.ba_tsteps           = [3, 8, 17, 21];             %Row vec of tsteps to use in BA
options.ba_pts              = [2, 3, 10, 16, 20, 23];     %Row vec of pts to use in BA

%Stereo Options
options.stereo.pts          = options.pts;
options.stereo.cams         = [301:305];
options.stereo.tstart       = 1;
options.stereo.tstop        = 10;

%Trajectory Estimation Options
% options.est.pts             = options.pts;
% options.est.tstart          = 1;
% options.est.tstop           = options.tstop - options.tstart+1;

%Plot Options
% options.plot.pts            = options.pts;
% options.plot.reprojframe    = 220;
options.plot.tstart         = 1;
options.plot.tstop          = (options.tstop - options.tstart)-(options.plot.tstart-1);
options.plot.linestyle1     = {'.-r','.-b','.-g', '.-m','.-k','.-c','.--r','.--b','.--g','+-r','+-b','+-g', '+-m','+-k','+-c','+--r','+--b','+--g'};
options.plot.linestyle2     = {'+-r','+-b','+-g', '+-m','+-k','+-c','+--r','+--b','+--g','o-r','o-b','o-g', 'o-m','o-k','o-c','o--r','o--b','o--g'};
options.plot.linestyle3     = {'o-r','o-b','o-g', 'o-m','o-k','o-c','o--r','o--b','o--g','.-r','.-b','.-g', '.-m','.-k','.-c','.--r','.--b','.--g'};
options.plot.colors         = {'r', 'g', 'b', 'c', 'm', 'k'};
%options.plot.savepath       = 'D:\Users\Matt\Documents\VT\Research\Motion_estimation_dev\Papers_and_Presentations\Bender2015scitech\';
options.plot.savefig        = 0;
options.plot.saveim_reproj  = 0;
options.plot.saveim_reproje = 0;
options.plot.fig_txt_props  = {'FontName', 'Times New Roman', 'FontSize', 18, 'FontWeight', 'Bold'};

%% Run the Import Utility 
Cam = load_svoboda_cal3(options);
