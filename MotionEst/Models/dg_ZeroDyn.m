function Gj_k = dg_ZeroDyn(x_km1,u_k,t)
%MMJACZERODYN       -Computes the jacobian of the motion model for zero
%dynamics given the last state, current controls, and current time step.

n = length(x_km1);
Gj_k = eye(n,n);

end