function derive_msmt_model()

syms f1 f2 a c1 c2 real;
syms om1 om2 om3 t1 t2 t3 real;
syms X Y Z real;

% IMPORTANT: if norm(om) = 0, need to set om such that it is a full 2*pi
% rotation without norm(om) = 0 inside calling code.
% (For example: om = [2*pi; 0; 0] or om = 2*pi/sqrt(3)*[1; 1; 1])
% IMPORTANT
om = [om1 om2 om3]';
theta = norm(om); 
omega = om/theta;
alpha = cos(theta);
beta = sin(theta);
gamma = 1-cos(theta);
omegav=[[0 -omega(3) omega(2)];
        [omega(3) 0 -omega(1)];
        [-omega(2) omega(1) 0 ]];
A = omega*omega';
R = eye(3)*alpha + omegav*beta + A*gamma;

t = [t1 t2 t3]';
Hinv = [  R   t;
        0 0 0 1];

p = [f1 f2 a c1 c2 om1 om2 om3 t1 t2 t3]';


K = [f1  a*f1   c1  0; 
     0    f2    c2  0;
     0    0      1  0;
     0    0      0  1];

M = 1/([0 0 1 0]*Hinv*[X; Y; Z; 1])*Hinv*[X; Y; Z; 1];

disp('deriving h(X)...')
h_sym = [1 0 0 0;
         0 1 0 0]*K*M;
% disp('deriving dh_dX...')
% dh_dx_sym = jacobian(h_sym,[X; Y; Z]);
disp('deriving dh_dp...')
dh_dp_sym = jacobian(h_sym,p);


% disp('writing h(X) to file...')
% h = matlabFunction(h_sym, 'vars', {[p; X; Y; Z]}, 'file', 'h_mfile');
% disp('writing dh_dX to file...')
% dh_dx = matlabFunction(dh_dx_sym, 'vars', {[p; X; Y; Z]}, 'file', 'dhdx_mfile');
disp('writing dh_dp to file...')
dh_dp = matlabFunction(dh_dp_sym, 'vars', {[p; X; Y; Z]}, 'file', 'dhdp_mfile');
disp('done')