function H = HTransform(q,link)
nstates = length(q);
%groups  = unique([link.Group]);
%g_max   = max(groups);
links = [];
nDof = 0;
for ll = 1:length(link)
    links = [links, ll];
    nDof = nDof + link(ll).nDof;
    if nDof==nstates
        break
    end
end
H = eye(4,4);

jj = 1;
path = kinpath(link,links(1),links(end));
q = q([link(path).StateInds],1);
for ll = path
    thetas = zeros(link(ll).nDof,1);
    disps  = zeros(link(ll).nDof,1);
    %if identification kernal is DH -> pull out vars accordingly
    if strcmp(link(ll).IDkern,'DH')
        if ll>1
            if strcmp(link(ll-1).IDkern,'YPR')
                H = if_last_ypr(H,ll,link);
            end
        end
        %for all degrees of freedom for this link
        for ii = 1:link(ll).nDof
            %if dof is translational 
            if link(ll).tDof(ii)
                thetas(ii,1)  = 0;
                disps(ii,1) = q(jj);
            %else dof is rotational
            else
                thetas(ii,1) = q(jj);
                disps(ii,1) = 0;
            end
            jj = jj+1;
        end
        %Form H based on DH convention and store in tree
        H = H*DHTransforms(link(ll).thetas+thetas,link(ll).alphas,...
                               link(ll).disps+disps,link(ll).offsets);
    %if identification kernal is YPR, -> pull out vars accordingly
    elseif  strcmp(link(ll).IDkern,'YPR') 

        trans = 1;
        rot   = 1;
        %for all degrees of freedom for this link
        for ii = 1:link(ll).nDof
            %if dof is translational
            if link(ll).tDof(ii)
                disps(trans,1) = q(jj);
                trans = trans+1;
            %if dof is rotational
            else
                thetas(rot,1) = q(jj);
                rot = rot+1;
            end
            jj = jj+1;
        end
        %Form H based on YPR convention and store in tree
        H = H*YPRTransform(link(ll).thetas+thetas, link(ll).disps+disps);
    end
end

function H = if_last_ypr(H,ll,link)

H = H+[zeros(4,3),[H(1:3,1:3)*link(ll-1).BFvecs(:,link(ll).ConPt);0]];


