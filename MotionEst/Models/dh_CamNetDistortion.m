function Hj_k = dh_CamNet(x_k, Ring)

npts = length(x_k)/3;
nrings = length(Ring);
nc_r = length(Ring(1).Cam);
Pi0 = [1,0,0,0;0,1,0,0];
Hj_k = zeros(2*npts*nc_r*nrings,3*npts);
z_hat = [0;0;1;0];
%Grab a mean
for rr = 1:nrings           %for each ring
    for cc = 1:nc_r         %for each camera
        for pp = 1:npts
            x_barkpp = [x_k(3*(pp-1)+1:3*pp);1];
            Hin = invH(Ring(rr).Cam(cc).H);
            %Determine predicted range to point
            lambda = z_hat'*Hin*x_barkpp; 
            %Determine Sensor Model Jacobian
            dlam_dx = z_hat'*Hin*(lambda)^(-2);
            %Determine partial of range wrt x
            Hj_k((rr-1)*nc_r*2*npts+(cc-1)*2*npts+2*(pp-1)+1:(rr-1)*nc_r*2*npts+(cc-1)*2*npts+2*pp,3*(pp-1)+1:3*pp) = (1/lambda*Pi0*[Ring(rr).Cam(cc).K,[0;0;0];0,0,0,1]*Hin-Pi0*[Ring(rr).Cam(cc).K,[0;0;0];0,0,0,1]*Hin*x_barkpp*dlam_dx)*[eye(3);0,0,0];
        end
    end
end