function options = create_meas_vec(options)
%This function creates the state vector indicies of the chain definition
%contianed within OPTIONS.LINK.   State vector indicies for the DOF
%variables of each link are stored in OPTIONS.LINK(ll).StateInd. Some links
%must be identified in groups, so the vector starts with group 1 and ends
%with group G.  

link = options.link;
nlink = length(link);

groups = [];
for ll = 1:nlink
    groups = [groups, link(ll).Group];
end
groups = unique(groups);

indx = 1;
for gg = groups
    for ll = get_group_links(link,gg);
        nvecs = size(link(ll).BFvecs,2);
        link(ll).MeasInds = indx:(indx+2*nvecs-1);
        indx = indx + 2*nvecs;
    end
end

options.nmeas = indx-1;
options.link = link;






    