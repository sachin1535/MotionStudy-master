function [y_bark, lambdas] = CamNetDistortion(x_k, camstruct)

cams = [];
for cc = 1:length(camstruct)
   if ~isempty(camstruct(cc).pts)
       cams = [cams, cc];
   end
end

[rows, ~] = size(x_k);
npts = rows/3;
ncam = length(cams);
%Pi0 = [1,0,0,0;0,1,0,0];
z_hat = [0;0;1;0];
%Grab a mean
y_bark = zeros(2*ncam*npts,1);
lambdas = zeros(npts,1,ncam);

    for cc = 1:ncam         %for each camera
        for pp = 1:npts
            dist_c = camstruct(cams(cc)).kc;
            x_barkpp = [x_k(3*(pp-1)+1:3*pp);1];
            Hin = invH(camstruct(cams(cc)).H);
            %Determine predicted range to point
            lambdas(pp,1,cc) = z_hat'*Hin*x_barkpp; 
            if lambdas(pp,1,cc) <0
                lambdas(pp,1,cc) = NaN;
            end
            x_cam = Hin*x_barkpp;
            x_n = x_cam(1:2)/lambdas(pp,1,cc);
            r = norm(x_n);
            
            x_d = (1+dist_c(1)*r^2+dist_c(2)*r^4+dist_c(5)*r^6)*x_n+[2*dist_c(3)*x_n(1)*x_n(2)+ dist_c(4)*(r^2+2*x_n(1)^2);2*dist_c(4)*x_n(1)*x_n(2)+dist_c(3)*(r^2+2*x_n(2)^2)];
            x_p = camstruct(cams(cc)).K*[x_d;1];
            y_bark((cc-1)*2*npts+2*(pp-1)+1:(cc-1)*2*npts+2*pp)=x_p(1:2);
        end
    end
