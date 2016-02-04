function [x_k, Rk] = ZeroDyn(x_km1,u_k,t, Rt_handle)
%MMJACZERODYN       -Computes the jacobian of the motion model for zero
%dynamics given the last state, current controls, and current time step.

n = length(x_km1);
x_k = x_km1;

Rt = Rt_handle();
Rk_cell = repmat({Rt},1,n/3);
Rk = blkdiag(Rk_cell{:});
% Rk = eye(n);

end