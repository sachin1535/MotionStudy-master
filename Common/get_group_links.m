function links = get_group_links(link, group)
%returns the link numbers in a vector LINKS from the structure LINK in
%group number GROUP

nlinks = length(link);
links = [];
for gg = group
    for ll = 1:nlinks
        if gg == link(ll).Group
            links = [links,ll];
        end
    end
end