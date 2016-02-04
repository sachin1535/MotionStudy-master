function plot_kin_chain(kinc,kinConfig, tsteps)
%PLOT_KIN_CHAIN     -plots the links of the kinematic chain defined by the
%structure CHAIN at the timesteps TSTEPS.  The kinematic chain will be
%plotted in the zero basis of link 1.
hold on
nlinks = length(kinConfig.link);
for ll = 1:nlinks
    points = [kinConfig.link(ll).BFvecs, kinConfig.link(ll).BFvecs(:,1)];
    n_bf_pts = size(points,2);
    %n_pts_tot = n_pts_tot+n_bf_pts;
    for tt = tsteps
        X = [eye(3,3),zeros(3,1)]*hnode2node(kinc(tt),kinConfig,1,ll)*[points;ones(1,n_bf_pts)];
        if tt==tsteps(1)
            linespec = '-m';
        else
            linespec = '-k';
        end
        feat_manip = 1000*YPRTransform([0,-15*pi/180,-5/180*pi],[.300,1.200,.400])*[0,1,0,0;0,0,1,0;1,0,0,0;0,0,0,1]'*[X;ones(1,size(X,2))];
        plot3(feat_manip(1,:)',feat_manip(2,:)', feat_manip(3,:)', linespec, 'LineWidth', 2)
    end
end