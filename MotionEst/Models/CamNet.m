function [y_bark, Qk] = CamNet(x_k, camstruct)

cams = [];
for cc = 1:length(camstruct)
   if ~isempty(camstruct(cc).H)
       cams = [cams, cc];
   end
end

[rows, cols] = size(x_k);
npts = rows/3;
ncam = length(cams);
Pi0 = [1,0,0,0;0,1,0,0];
z_hat = [0;0;1;0];

%Grab a mean
y_bark = zeros(2*ncam*npts,1);
Qk = eye(length(y_bark));
for cc = 1:ncam         %for each camera
    for pp = 1:npts
        x_barkpp = [x_k(3*(pp-1)+1:3*pp);1];
        Hin = invH(camstruct(cams(cc)).H);
        %Determine predicted range to point
        lambda = z_hat'*Hin*x_barkpp; 
        %Determine Sensor Model Jacobian
        ndx = (cc-1)*2*npts+(2*(pp-1)+1:2*pp);
        y_bark(ndx) = 1/lambda*Pi0*[camstruct(cams(cc)).K,[0;0;0];0,0,0,1]*Hin*x_barkpp;
%         p = [camstruct(cams(cc)).foc_l;
%              camstruct(cams(cc)).skew; 
%              camstruct(cams(cc)).prin_p;
%              camstruct(cams(cc)).om;
%              camstruct(cams(cc)).T];
%         sigpvec = [(camstruct(cams(cc)).foc_l_e/3).^2;
%                    (camstruct(cams(cc)).skew_e/3).^2; 
%                    (camstruct(cams(cc)).prin_p_e/3).^2;
%                    (camstruct(cams(cc)).om_e/3).^2;
%                    (camstruct(cams(cc)).T_e/3).^2];
%         Jp = dhdp_mfile([p; x_barkpp]);
%            
%         Q2 = Jp*diag(sigpvec)*Jp';
% 
%         uncertainty_increase_factor = 1.2;
%         % multiply by this to increase uncertainty due to unmodeled things
%         % (uncertainty in distortion, higher-order-uncertainty, etc)
%         Qk(ndx,ndx) = Q2*uncertainty_increase_factor;
    end
end
