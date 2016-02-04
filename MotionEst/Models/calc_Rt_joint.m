function Rt = calc_Rt_joint(ll, link)
% This function generates the additive uncertainty for the motion model in
% joint space.  The output is a [nDof X nDof] covariance matrix where nDof
% is the number of degrees of freedom of the link LL for which Rt is
% calculated.  If LL is a vector of link numbers, all links are included in
% the covarience matric

%units are radians for angle and mm for distance
nDof = sum([link(ll).nDof]);
tDof = [link(ll).tDof];
Rt = zeros(nDof,nDof);
%%%%Future, consider increasing uncertainty inteligently for groups with
%%%%multiple links*******************
for dof = 1:nDof
    if tDof(dof)
        Rt(dof,dof) = .20^2;           %assume position uncertainty of 10mm
    else
        Rt(dof,dof) = (7*pi/180)^2;     %assume position uncertainty equivalent to 5 degrees
    end
end
        