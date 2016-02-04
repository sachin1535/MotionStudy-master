function nodes = kinpath(tree, n1, n2)
%KINPATH        -Creates a row vector of nodes numbers required to travel
%from node 1 to node 2.  Inputs are a kinematic tree structure (tree), a
%start node (n1) and a stop node (n2).

%save start and stop nodes to a dynamic node variable
node1 = n1;
node2 = n2;
%seed a back appended node array for traveling from the start node
nodesf = n1;
%seed a front appended node array for traveling from the stop node
nodesb = n2;

while node1~=node2                      %while node1 and node2 are not equal
    if node1<node2                      %if node1 is less than node2
        node2 = tree(node2).parent;     %go to parent of node2
        nodesb = [node2, nodesb];       %append new node2 in front of current
    else                                %else, node2 is less than node1
        node1 = tree(node1).parent;     %go to parent of node1
        nodesf = [nodesf, node1];       %append new node1 behind current 
    end
end
l1 = length(nodesf);                    %check length of each array
l2 = length(nodesb);
if l1>l2                                %chop off one node number from longest array. 
    nodesf = nodesf(1:l1-1);
else
    nodesb = nodesb(2:l2);
end
nodes = [nodesf, nodesb];               %Concatenate the arrays