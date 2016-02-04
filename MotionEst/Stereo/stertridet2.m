function x = stertridet2(phi1, phi2, camstruct1,camstruct2)
K1 = camstruct1.K;
K2 = camstruct2.K;
H1 = camstruct1.H;
H2 = camstruct2.H;

%Determine Vector from o1 to o2 in world frame.
P1= H1(1:3,4);
P2 = H2(1:3,4);
T = P2-P1;

%determine unit vectors to observed by cameras in world frame
u1 = H1(1:3,1:3)*(K1\[phi1;1])/norm((K1\[phi1;1]));
u2 = H2(1:3,1:3)*(K2\[phi2;1])/norm((K2\[phi2;1]));

%Determine Plane parallel to both lines. 
N_hat = cross(u1,u2);

%Project vector between coordinate frames onto normal vector (i.e. the distance between planes)
d_vec = T'*N_hat*N_hat;
%dotN_hat_T = 
%d_vec = d;

T_til = T-d_vec;

gamma = acos(u1'*u2/(norm(u1)*norm(u2)));
beta  = pi-acos(T_til'*u2/(norm(T_til)*norm(u2)));
alpha = acos(T_til'*u1/(norm(T_til)*norm(u1)));
b = norm(T_til)/sin(gamma)*sin(beta)*u1+P1;
a = norm(T_til)/sin(gamma)*sin(alpha)*u2+P2;

%d_vec_new = b-(a+T);
x = b/2+a/2;

% plot3(x(1),x(2),x(3),'ob')
% plot3(a(1),a(2),a(3),'og')
% plot3(b(1),b(2),b(3),'or')

