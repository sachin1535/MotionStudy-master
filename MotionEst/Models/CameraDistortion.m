function xp = CameraDistortion(x, camstruct)

cams = [];
for cc = 1:length(camstruct)
   if ~isempty(camstruct(cc).pts)
       cams = [cams, cc];
   end
end

ncam = length(cams);
npts = length(x)/(2*ncam);

    for cc = 1:ncam         %for each camera
        K = camstruct(cams(cc)).K;
        dist_c = camstruct(cams(cc)).kc;
        for pp = 1:npts
            x_n = [1 0 0 ; 0 1 0]*inv(K)*[x((cc-1)*2*npts+2*(pp-1)+1:(cc-1)*2*npts+2*pp);1];
            r = norm(x_n);

            xd = (1+dist_c(1)*r^2+dist_c(2)*r^4+dist_c(5)*r^6)*x_n+[2*dist_c(3)*x_n(1)*x_n(2)+ dist_c(4)*(r^2+2*x_n(1)^2);2*dist_c(4)*x_n(1)*x_n(2)+dist_c(3)*(r^2+2*x_n(2)^2)];

            xp((cc-1)*2*npts+2*(pp-1)+1:(cc-1)*2*npts+2*pp) = K(1:2,:)*[xd;1];
        end
    end
end