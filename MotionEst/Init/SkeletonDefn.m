%%BROWNKINDEF12         -This program was created to solve the motion
%%identification problem for flapping wing flight in bats.  While the code
%%below was specifically created for this problem, it can be adapted to
%%solve any motion identification problem for an open kinematic chain.  The
%%data supplied to this program is 3D point trajectories of all joint links
%%and robotic configuration information in the form of a DH table.  The
%%output is graphs of the joint variable for each degree of freedom, plots 
%%of motion identificaiton verification (i.e. experimental vs. identified
%%trajectories in cartesian coordinates, and a plot of all experimental
%%points and all identified points on the kinematic chain for each time
%%step.

%%Future revisions of this code will contain the following:
%%------1. GUI based kinematic configuration definition.  Currently this is
%%done manually, but a GUI will make the process easier and less prone to
%%bugs. 

%%------2. GUI based data import.  Data will be associated to the correct
%%links and a configuration file will be created for future processing so
%%that the definition does not have to be performed for each run.

%%------3. Trajectory smoothing.  All of the identified trajectories 

%% Define Kinematic Chain Structure For Motion Identification Software
%This program is a manual assignment of the hierarchy implemented for the
%motion identification code which will follow.  It assigns node numbers and
%string identifiers, as well as the parent child relationships between
%nodes.  In the future this portion of the motion identification software
%will be automated.  


% Node Names and Parent-Child Info
%             Node Name   Parent    Child       Connecting Point       Group   ID Kernal  
NodeNames   = {'BB'      ,   [] ,    [2],            [],                1,      'YPR';
              'RHum'     ,   [1],    [3],            [3],               2,      'DH';
              'RRad'     ,   [2],    [4,5,6],        [2],               2,      'DH';
              'RD3Met'   ,   [3],    [],             [2],               3,      'DH'; 
              'RD4Met'   ,   [3],    [],             [2],               3,      'DH';
              'RD5Met'   ,   [3],    [],             [2],               3,      'DH'};
          
nnodes = length(NodeNames);

for nn = 1:nnodes 
    
synthConfig.link(nn).nnames = NodeNames{nn,1};
synthConfig.link(nn).parent = NodeNames{nn,2};
synthConfig.link(nn).child  = NodeNames{nn,3};
synthConfig.link(nn).ConPt  = NodeNames{nn,4};
synthConfig.link(nn).Group  = NodeNames{nn,5};
synthConfig.link(nn).IDkern = NodeNames{nn,6};

end


%% Create DOF Tree
DOFAssign = {'BB', 6; 'Hum', 3; 'Rad', 1; 'Met', 2; 'Phal', 1};

%Cycle through labels and assign NDOFs
for ii = 1:length(DOFAssign);
    %Get current node label 
    cmpstr = strcat('.*',DOFAssign(ii,1),'.*');
    nlist = [];
    %find nodes with label identifer
    for kk = 1:nnodes
        TF = regexp(synthConfig.link(kk).nnames, cmpstr{1});
        %if node label contains identifier
        if TF
            %add node to list
            nlist = [nlist, kk];
        end
    end
    %for all matched nodes
    for jj = nlist
        %assign number of Dofs
        synthConfig.link(jj).nDof = DOFAssign{ii,2};
    end
    
end

%Dof Type Specification (1 for translational, 0 for rotational)
for ii = 1:length(synthConfig.link);
    if synthConfig.link(ii).nDof == 6;
        synthConfig.link(ii).tDof = [1;1;1;0;0;0]';
    else
        synthConfig.link(ii).tDof = zeros(synthConfig.link(ii).nDof,1)';
    end
end


%% Define DH params
nn = 1;
%----------------------------------Base Body CF Defn---------------------------------------
synthConfig.link(nn).thetas  = [0;0;0;0;0;0];
synthConfig.link(nn).alphas  = [pi/2; pi/2; 0; pi/2; pi/2; 0];
synthConfig.link(nn).disps   = [0;0;0;0;0;0];
synthConfig.link(nn).offsets = [0;0;0;0;0;0];
synthConfig.link(nn).H = DHTransforms(synthConfig.link(nn).thetas,synthConfig.link(nn).alphas,synthConfig.link(nn).disps,synthConfig.link(nn).offsets);
nn = nn+1;
%----------------------------------Humerus CF Defn---------------------------------------
hum_len = 60;
synthConfig.link(nn).thetas  = [0;0;0];
synthConfig.link(nn).alphas  = [pi/2; -pi/2; pi/2];
synthConfig.link(nn).disps   = [0;0;hum_len];
synthConfig.link(nn).offsets = [0;0;0];
synthConfig.link(nn).H = DHTransforms(synthConfig.link(nn).thetas,synthConfig.link(nn).alphas,synthConfig.link(nn).disps,synthConfig.link(nn).offsets);
nn = nn+1;
%----------------------------------Raduis CF Defn---------------------------------------
gamma = 85/180*pi;
synthConfig.link(nn).thetas  = [0];
synthConfig.link(nn).alphas  = [0];
synthConfig.link(nn).disps   = [0];
synthConfig.link(nn).offsets = 70;
synthConfig.link(nn).H       = DHTransforms(synthConfig.link(nn).thetas,synthConfig.link(nn).alphas,synthConfig.link(nn).disps,synthConfig.link(nn).offsets);
nn = nn+1;

%----------------------------------Digit 3 Metacarpal CF Defn---------------------------------------
%met3_len   = LinkLen(browne,nn);
synthConfig.link(nn).thetas  = [0;0];
synthConfig.link(nn).alphas  = [pi/2;0];
synthConfig.link(nn).disps   = [0;0];
synthConfig.link(nn).offsets = [0;70];
synthConfig.link(nn).H = DHTransforms(synthConfig.link(nn).thetas,synthConfig.link(nn).alphas,synthConfig.link(nn).disps,synthConfig.link(nn).offsets);
nn = nn+1;

%----------------------------------Digit 4 Metacarpal CF Defn---------------------------------------
%met4_len   = LinkLen(browne,nn);
synthConfig.link(nn).thetas  = [0;0];
synthConfig.link(nn).alphas  = [pi/2;0];
synthConfig.link(nn).disps   = [0;0];
synthConfig.link(nn).offsets = [0;70];
synthConfig.link(nn).H = DHTransforms(synthConfig.link(nn).thetas,synthConfig.link(nn).alphas,synthConfig.link(nn).disps,synthConfig.link(nn).offsets);
nn = nn+1;
%----------------------------------Digit 5 Metacarpal CF Defn---------------------------------------
%met5_len   = LinkLen(browne,nn);
synthConfig.link(nn).thetas  = [0;0];
synthConfig.link(nn).alphas  = [pi/2;0];
synthConfig.link(nn).disps   = [0;0];
synthConfig.link(nn).offsets = [0;70];
synthConfig.link(nn).H = DHTransforms(synthConfig.link(nn).thetas,synthConfig.link(nn).alphas,synthConfig.link(nn).disps,synthConfig.link(nn).offsets);
nn = nn+1;

%% Create an articulated mechanism
%Create Body Fixed Vectors, parent child relations, and DH params
%link 1 
synthConfig.link(1).BFvecs(:,1) = [0,0,0]';
synthConfig.link(1).BFvecs(:,2) = [35,10,0];
synthConfig.link(1).BFvecs(:,3) = [35,-10,0];
synthConfig.link(1).BFvecs(:,4) = [30,0,-5];
synthConfig.link(1).BFvecs(:,5) = [30,0,5];

%link 2
synthConfig.link(2).BFvecs(:,1) = [0,0,0]';
synthConfig.link(2).BFvecs(:,2) = [40,0,0];

%link 3
synthConfig.link(3).BFvecs(:,1) = [0,0,0]';
synthConfig.link(3).BFvecs(:,2) = [60,0,0];

%link 4
synthConfig.link(4).BFvecs(:,1) = [0,0,0]';
synthConfig.link(4).BFvecs(:,2) = [60,0,0];

%link 5
synthConfig.link(5).BFvecs(:,1) = [0,0,0]';
synthConfig.link(5).BFvecs(:,2) = [60,0,0];

%link 6
synthConfig.link(6).BFvecs(:,1) = [0,0,0]';
synthConfig.link(6).BFvecs(:,2) = [60,0,0];







