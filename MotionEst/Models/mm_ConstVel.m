function [x_k, Rk] = mm_ConstVel(x_km1, x_km2, options, links, Rt_handle)
%MMJACZERODYN       -Computes the jacobian of the motion model for zero
%dynamics given the last state, current controls, and current time step.
link = options.link;
%links = get_group_links(link,group);
%l_max = max(links);
%ind_max = max(link(l_max).StateInds);
%num_states_comp = length(x_k);
%A1 = eye(num_states_comp,num_states_comp);
% A2 = zeros(ind_max,length(x_km2));
% x = zeros(ind_max+options.nstate,1);
% %x(1:num_states_comp) = x_k;
% %A1 = 2*eye(length([link(links).StateInds]),length([link(links).StateInds]);
% for ll = links
%     %A1(link(ll).StateInds-link(links(1)).StateInds(1)+1,link(ll).StateInds-link(links(1)).StateInds(1)+1) = 2*eye(length([link(ll).StateInds]));
%     A2(link(ll).StateInds,options.nstate+link(ll).StateInds) = -1*eye(length(link(ll).StateInds),length(link(ll).StateInds));
%     x(link(ll).StateInds) = x_km1(link(ll).StateInds);
% end
% x(ind_max+1:ind_max+1+options.nstate) = x_km2;

%x_k = 2*x_km1 - x_km2;
x_k = x_km1;
Rk = Rt_handle(links);
%Rk = [zeros(num_states_comp,ind_max);
%      zeros(ind_max-num_states_comp,ind_max-num_states_comp),Rk];

end