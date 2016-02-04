function [options, camstruct] = load_svoboda_cal3(camstruct, options)
%determine the cameras to be imported
[options.ext_cal_file,options.ext_cal_path]  = uigetfile;
load([options.ext_cal_path,options.ext_cal_file],'C','R','cam');
fprintf('Extrinsic Parameters are Available for the Following Cameras:\n')
for cc = 1:length(cam)
    fprintf('%d\t',cam(cc).camId);
end
cams = input('\n Which cams would you like to import extrinsic parameters for?:');
%determine the points to be estimated

%default_dir = pwd;
%plot_flag = options.plotflag;

%Determine the files in this directory
%list = dir(options.path);
%initialize index vars
%mm = 1;

%if the global frame is 1 then the 

for cc = cams
    index = find([cam.camId]==cc);
    cam_num = cc;
    fid = fopen([options.ext_cal_path,filesep,'Cam',num2str(cam_num),'.cal']);
    data = textscan(fid,'%s%s%f');
    %t = data{3}(10:12);
    t = -C(:,index);
    %R = [data{3}(1:3)';data{3}(4:6)';data{3}(7:9)'];
    Rot = R(3*index-2:3*index,:)';
    camstruct(cam_num).H = [Rot,t;0,0,0,1];
    camstruct(cam_num).K = [data{3}(13:15)';data{3}(16:18)';data{3}(19:21)'];
    for ii = 1:9
        if camstruct(cam_num).K(ii)<0
            camstruct(cam_num).K(ii) = -camstruct(cam_num).K(ii);
        end
    end

    fclose('all');
end
