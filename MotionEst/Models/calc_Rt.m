function Rt = calc_Rt()
% This function generates the additive uncertainty for the motion model
% x_t = eye()*x_{t-1} for bat-flight points based on some
% back-of-the-envelope assumptions about how far each point could move
% per frame.  The idea is to construct a 95% confidence ellipse for the
% new location of the point based on a rough approximation of it's
% maximum movement in a single frame.

% A number of assumptions are inherent in this model which may be
% revisited or revised, I hope those are evident from the names.  In
% brief, this model assumes rigid-wing flapping and calculates the
% displacement of a point on the tip of the rigid wing, then square-sums
% that with the displacement from the bat-flight through the world.

wingspan = .5; %[m] (from wingtip to wingtip, divide by 2 for single wing)
flap_freq = 4; %[Hz]
framerate = 120; %[fps]
flap_angle = pi/2; %[rad] (only downstroke, mult by 2 for total down-and-up-flap)
bs_factor = 1;

flapdist_per_frame = flap_angle*2*(wingspan/2)*(flap_freq)*(1/framerate)*bs_factor;
flapdist_per_frame = flapdist_per_frame * 1000; % convert to [mm]

flight_speed = 3; %[m/s]
flydist_per_frame = flight_speed*(1/framerate);
flydist_per_frame = flydist_per_frame * 1000; % convert to [mm]

uncertainty_radius = sqrt(flapdist_per_frame^2 + flydist_per_frame^2);
%uncertainty_radius = 50;
% chi2inv is from MATLAB statistics toolbox.  In case you don't have it,
% here are a few useful values, just write them in instead of chi2inv
% function:
% chi2inv(.90,3) = 6.2514
% chi2inv(.95,3) = 7.8147
% chi2inv(.98,3) = 9.8374
% chi2inv(.99,3) = 11.3449
uncertainty_ratio = .95;
% variance = (uncertainty_radius^2)/chi2inv(uncertainty_ratio,3);
variance = (uncertainty_radius^2)/7.8147;
Rt = variance*eye(3);
%Rt = eye(size(Rt));