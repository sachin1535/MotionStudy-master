function camstruct = load_caltech_intrinsic(camstruct, options, cams)
%LOAD_CALTECH_INTRINSIC     -loads the intrinsic parameters from the
%CalTech calibration toolbox.  Parameters must be stored in a MAT file
%which contains KK (intrinsic camera parameters) and KC (vector of
%distortion coefficients.

%set the location to get the parameters from 
path = [options.path,filesep,'..',filesep,'Calibration_run',filesep,'Intrinsic',filesep,'CalTech'];
%for the desired cameras, import the parameters
for cc = cams
    load([path,filesep,'Cam',num2str(cc),filesep,'int_cam',num2str(cc),'.mat'],'KK','kc')
    camstruct(cc).K = KK;
    camstruct(cc).kc = kc;
end

    