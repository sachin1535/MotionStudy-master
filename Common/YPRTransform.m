function [ H ] = YPRTransform(theta,d)
%% funrollpitchyaw :  calculate the rotation matrix for 3-2-1 Euler angles
%
%   Inputs :    
%   roll    -   roll angle in radians, rotation about 1 axis
%   pitch   -   pitch angle in radians, rotation about 2 axis
%   yaw     -   yaw angle in radians, rotation about 3 axis
%
%%   yaw about  3 axis
%

roll  = theta(1);
pitch = theta(2);
yaw   = theta(3);

dx = d(1);
dy = d(2);
dz = d(3);

R3=[cos(yaw),-sin(yaw),0;
    sin(yaw),cos(yaw),0;
    0       ,0       ,1];

%
%%   pitch about 2 axis
%
R2=[cos(pitch)  ,       0,  sin(pitch);
    0           ,       1,  0;
    -sin(pitch) ,       0,  cos(pitch)];
%
%%   roll about 1 axis
%

R1=[1,  0,           0;
    0,  cos(roll),  -sin(roll);
    0,  sin(roll),  cos(roll)];

%
%%   create the 3-2-1 rotation matrix
%

R=R3*R2*R1;
H = [R, [dx;dy;dz];0,0,0,1];

end

