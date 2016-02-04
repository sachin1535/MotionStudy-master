function H = hnode2node(kinD, kinConfig, node1, node2)
%HN2N       Creates homogeneous transform from the basis at node 1 to the 
%basis at node 2.  KINTREE is the kinematic tree structure defined for this 
%problem. NODE1 and NODE2 are string variables which specify the points to 
%move between. H is a single rotation matrix which expresses the body fixed 
%frame at node 1 in the body fixed frame at node 2.

%%Determine type of node identifier
if ischar(node1)    %if node 1 is a character string, find on map
    n1 = chainsearch(kinConfig.link, node1);
else                %else assume numeric value of node was entered
    n1 = node1;
end

if ischar(node2)    %if node 2 is a character string, find on map
    n2 = chainsearch(kinConfig.link, node2);
else                %else assume numerica value of node was entered
    n2 = node2; 
end

%%If node1 == node2, return H matrix at that node and exit
if n1==n2
    H = kinD.link(n1).H;
    return 
end

%determine path between nodes
path = kinpath(kinConfig.link,n1,n2);
if length(path)>1
    dir  = sign(path(2)-path(1));
else
    dir = 1;
end
minnode = min(path);
if minnode == path(1);
    sw = 0;
elseif minnode == path(length(path));
    sw = 0;
else
    sw = 1;
end

%Determine number of points in time series
ndata = length(kinD);

%Seed H output matrix with identity matricies 
H = zeros(4,4,ndata);
for ii = 1:ndata
H(:,:,ii) = eye(4,4);
end
node_last = path(1);
for ii = 1:ndata  %for each time step
    jj = 2;         
    for node = path %for all nodes on path
        node_next = path(jj);
        if dir <0  %if traveling from leaves to trunk, pre mult by H transform
            if ~strcmp(kinConfig.link(node_last).IDkern,kinConfig.link(node).IDkern)
                H(:,:,ii) = [H(1:3,1:3),H(1:3,1:3)*kinConfig.link(node).BFvecs(:,kinConfig.link(node_next).ConPt)+H(1:3,4);0,0,0,1];
            else
                H(:,:,ii) = invH(kinD(ii).link(node).H)*H(:,:,ii);
            end
        elseif dir > 0  %if traveling from trunk to leaves, pre-mult by transpose of H Transform
            if ~strcmp(kinConfig.link(node_next).IDkern,kinConfig.link(node).IDkern)
                H(:,:,ii) = H(:,:,ii)*(kinD.link(node).H+[zeros(3,3),kinD.link(node).H(1:3,1:3)*kinConfig.link(node).BFvecs(:,kinConfig.link(node_next).ConPt);0,0,0,0]);
            else
                H(:,:,ii) = H(:,:,ii)*kinD(ii).link(node).H; 
            end
        end
        %If there is a direction switch in the path and the node == the
        %switching node, multiply by the transpose of the same H Transform
        %***This assumes that direction change MUST happen at the lowest
        %node in the series, and the first direction must have been from
        %leaves to trunk
        if sw && node == minnode
            H(:,:,ii) = invH(Hall(:,:,ii))*H(:,:,ii); 
        end
        node_last = node;
        jj = jj+1;
        if jj>length(path)
            jj = jj-1;
        end
    end
end
    
    

