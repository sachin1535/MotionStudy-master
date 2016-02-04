function options = create_state_vec(options)
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
        nDof = link(ll).nDof;
        link(ll).StateInds = indx:(indx+nDof-1);
        indx = indx + nDof;
    end
end
options.nstate = indx-1;
options.link = link;



function links = get_group_links(link, group)
%returns the link numbers in a vector LINKS from the structure LINK in
%group number GROUP

nlinks = length(link);
links = [];
for ll = 1:nlinks
    if group == link(ll).Group
        links = [links,ll];
    end
end


    