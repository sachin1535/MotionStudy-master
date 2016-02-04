function Hj_k = dh_CamNet(x_k, camstruct)

cams = [];
for cc = 1:length(camstruct)
   if ~isempty(camstruct(cc).H)
       cams = [cams, cc];
   end
end

npts = length(x_k)/3;
ncam = length(cams);
Pi0 = [1,0,0,0;0,1,0,0];
Hj_k = zeros(2*npts*ncam,3*npts);
z_hat = [0;0;1;0];
%Grab a mean
for cc = 1:ncam         %for each camera
    for pp = 1:npts
        x_barkpp = [x_k(3*(pp-1)+1:3*pp);1];
        Hin = invH(camstruct(cams(cc)).H);
        %Determine predicted range to point
        lambda = z_hat'*Hin*x_barkpp; 
        %Determine Sensor Model Jacobian
        dlam_dx = z_hat'*Hin*(lambda)^(-2);
        %Determine partial of range wrt x
        rndx = (cc-1)*2*npts+2*(pp-1)+1:(cc-1)*2*npts+2*pp;
        cndx = 3*(pp-1)+1:3*pp;
        Hj_k(rndx,cndx) = (1/lambda*Pi0*[camstruct(cams(cc)).K,[0;0;0];0,0,0,1]*Hin-Pi0*[camstruct(cams(cc)).K,[0;0;0];0,0,0,1]*Hin*x_barkpp*dlam_dx)*[eye(3);0,0,0];
    end
end
