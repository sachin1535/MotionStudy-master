function camstruct_adjusted = bundle_adjust_pixel2(camstruct, options)
%BUNDLE_ADJUST      -This function can be run after an import utility
%which creates the Ring(1).Cam structure.  This code adjusts the external
%camera calibration by minimizing the error between triangulation performed
%with any pair of cameras.  

%If bundle adjustment is not desired, exit this function
if ~options.ba
    camstruct_adjusted = camstruct;
    return
end
cams = input('Which cams to compute or load bundle adjustment for?:');
ncam = length(cams);

%Save Current Path and Add to Working Dir
%default_dir = pwd;
%addpath(default_dir);
%Move to Location of Data
%cd(options.path);

%Compile Original Extrinsics and Plot
H = zeros(4,4,ncam);
for cc = 1:ncam
    H(:,:,cc) = camstruct(cams(cc)).H;
end
figure
hold on
CFPlot(H, .1)
axis equal

savefile = 'bundle_adjust.mat';
%Determine if Bundle Adjust Already extists
if (options.ba == 1 || options.ba ==3)
    [bundle_file, bundle_dir] = uigetfile;
    load([bundle_dir,filesep,bundle_file])
    camstruct_adjusted = camstruct;
    for cc = cams
        camstruct_adjusted(cc).H = camstruct_exp(cc).H;
        camstruct_adjusted(cc).T = camstruct_exp(cc).T;
        camstruct_adjusted(cc).om = camstruct_exp(cc).om;
    end
    clear camstruct_exp
    for cc = 1:ncam
        H(:,:,cc) = camstruct_adjusted(cams(cc)).H;
    end
    
    CFPlot(H, .1)
    title('Extrinics Before and After Loading Previous Bundle Adjustment')
    axis equal
    %cd(default_dir)
    
    if options.ba ==1
        return
    else
        camstruct = camstruct_adjusted;
    end
else
    camstruct_adjusted = camstruct;
end


% %Set options of fminsearch to plot function evaluations.
% options_fmin = optimset('PlotFcns', @optimplotfunccount);
% %find the parameters that minimize error
% delta_vec = fminsearch(@(delta_vec) fcost(camstruct, options, delta_vec), zeros(3*2*ncam,1), options_fmin);
% 
% %Reset extrinsics to new values
% nparam = length(delta_vec);
% delta_vec_rot = delta_vec(1:nparam/2);
% delta_vec_pos = delta_vec(nparam/2+1:nparam);
% %delta_vec_rot = delta_vec;
% %delta_vec_pos = delta_vec;
% 
% camstruct_adjusted = camstruct;
% for cc = 1:ncam
%     camstruct_adjusted(cams(cc)).H = camstruct_adjusted(cams(cc)).H*[rodrigues(delta_vec_rot(3*(cc-1)+1:3*cc,1)),delta_vec_pos(3*(cc-1)+1:3*cc,1);zeros(1,3),1];
%     %camstruct_adjusted(cams(cc)).H = camstruct_adjusted(cams(cc)).H*[eye(3),delta_vec_pos(3*(cc-1)+1:3*cc,1);zeros(1,3),1];
%     Hin = invH(camstruct_adjusted(cams(cc)).H);
%     camstruct_adjusted(cams(cc)).T = Hin(1:3,4);
%     camstruct_adjusted(cams(cc)).om = rodrigues(Hin(1:3,1:3));
%     
% end
% camstruct_exp = camstruct_adjusted;
% save(savefile, 'camstruct_exp')
% cd(default_dir);
% 
% %Compile New Extrisics and Show on Same Plot
% for cc = 1:ncam
%     H(:,:,cc) = camstruct_adjusted(cams(cc)).H;
% end
% 
% CFPlot(H, 50)
% title('Extrinics Before and After Adjustement')
% axis equal


function cost = fcost(camstruct, options, delta_vec)
%Grab Required Options

npts = length(options.ba_pts);
pts = options.ba_pts;
timesteps = options.ba_tsteps;
nsteps = length(options.ba_tsteps);
nparam = length(delta_vec);
cams = options.cams;
ncam = length(cams);

delta_vec_rot = delta_vec(1:nparam/2);
delta_vec_pos = delta_vec(nparam/2+1:nparam);
%delta_vec_rot = delta_vec;
%delta_vec_pos = delta_vec;

%Modify extrinsics of camera setup
for cc = 1:ncam
    %camstruct(cams(cc)).H = camstruct(cams(cc)).H*[eye(3),delta_vec_pos(3*(cc-1)+1:3*cc,1);zeros(1,3),1];
    camstruct(cams(cc)).H = camstruct(cams(cc)).H*[rodrigues(delta_vec_rot(3*(cc-1)+1:3*cc,1)),delta_vec_pos(3*(cc-1)+1:3*cc,1);zeros(1,3),1];
    Hin = invH(camstruct(cams(cc)).H);
    camstruct(cams(cc)).T = Hin(1:3,4);
    camstruct(cams(cc)).om = rodrigues(Hin(1:3,1:3));
    
    if norm(camstruct(cams(cc)).om) == 0;
        camstruct(cams(cc)).om = [2*pi;0;0];
    end
    
end

%determine pairs of cameras
npair = 0;
pair_list = [];
for ii = 1:ncam-1
    npair = npair + ii;
    for jj = ii+1:ncam
        pair_list = [pair_list; cams(ii),cams(jj)];
    end
end

%Determine stereo reconstructions of each pair of cams with new extrinsics
X_ster = zeros(3*npts, nsteps, npair);
phi_ster = zeros(2*npts*ncam,nsteps,npair);

for pair = 1:npair
    for kk = 1:nsteps
        for pp = 1:npts
            X_ster(3*(pp-1)+1:3*pp,kk,pair) = stertridet2(camstruct((pair_list(pair,1))).pts(:,timesteps(kk),pts(pp)), ...
                                               camstruct((pair_list(pair,2))).pts(:,timesteps(kk),pts(pp)), ...
                                               camstruct((pair_list(pair,1))),camstruct((pair_list(pair,2))));
            [phi_ster(:,kk,pair),~]         = CamNet(X_ster(:,kk,pair), camstruct);
            
        end
    end
end

phi_meas = zeros(2*npts*ncam,nsteps);
for cc = 1:ncam
    for pp = 1:npts
        for kk = 1:nsteps
            phi_meas(2*(cc-1)*npts+2*(pp-1)+1:2*(cc-1)*npts+2*pp,kk) = camstruct(cams(cc)).pts(:,timesteps(kk),pts(pp));
        end
    end
end

phi_error = zeros(2*npts*ncam,nsteps,npair);
for pair = 1:npair
    phi_error(:,:,pair) = (phi_ster(:,:,pair)-phi_meas).*(phi_ster(:,:,pair)-phi_meas);
end

occlusions = isnan(phi_error);
phi_error(occlusions) = 0;

if any(isnan(phi_error))
    keyboard
end
w1 = 0;
w2 = 0;
%cost = sum(sum(sum(phi_error)))+w*(delta_vec'*delta_vec);
%cost = sum(sum(sum(phi_error)))+w2*(delta_vec_pos'*delta_vec_pos);
cost = sum(sum(sum(phi_error)))+w1*(delta_vec_rot'*delta_vec_rot)+w2*(delta_vec_pos'*delta_vec_pos);



