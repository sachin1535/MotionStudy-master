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

fprintf('Creating Kinematic Definition...\n')
% Node Names and Parent-Child Info
%             Node Name   Parent    Child       Connecting Point       Group   ID Kernal  
NodeNames   = {'BB'      ,   [] ,    [],            [],                1,      'YPR';
               'LHum'     ,   [1],    [3],            [3],               2,      'DH';
               'LRad'     ,   [2],    [4,5,6],        [2],               2,      'DH';
               'RD3Met'   ,   [3],    [],             [2],               3,      'DH'; 
               'RD4Met'   ,   [3],    [],             [2],               3,      'DH';
               'RD5Met'   ,   [3],    [],             [2],               3,      'DH'};
          
nnodes = size(NodeNames,1);

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

fprintf('---------------------------------------------------\n')
fprintf(' Link \t Name \t nDof \t Group \t Parent \tIdKern \n')
fprintf('---------------------------------------------------\n')

for ii = 1:length(synthConfig.link)
    fprintf(' %d \t\t %s',ii, synthConfig.link(ii).nnames)
    if length(synthConfig.link(ii).nnames)<3
        fprintf('\t\t %d \t\t\t %d \t\t %d \t\t %s \n', synthConfig.link(ii).nDof,synthConfig.link(ii).Group,synthConfig.link(ii).parent,synthConfig.link(ii).IDkern)
    elseif length(synthConfig.link(ii).nnames)<6
        fprintf('\t %d \t\t\t %d \t\t %d \t\t %s \n', synthConfig.link(ii).nDof,synthConfig.link(ii).Group,synthConfig.link(ii).parent,synthConfig.link(ii).IDkern)
    else
        fprintf('  %d \t\t\t %d \t\t %d \t\t %s \n', synthConfig.link(ii).nDof,synthConfig.link(ii).Group,synthConfig.link(ii).parent,synthConfig.link(ii).IDkern)
    end
end
fprintf('---------------------------------------------------\n')


%% Create Body Fixed Vectors
%specifiy the points which are on each link
synthConfig.link(1).pt_nums = [1,4,5];
synthConfig.link(2).pt_nums = [5,6];
synthConfig.link(3).pt_nums = [6,7,8];
synthConfig.link(4).pt_nums = [8,9,10];
synthConfig.link(5).pt_nums = [8,13,14];
synthConfig.link(6).pt_nums = [8,17];
%load the stereo triangulation data
load([options.path,filesep,'StereoStruct.mat']);
npair = length(Stereo);         %determine the number of cameras pairs

for ll = 1:length(synthConfig.link)
    vectors = [];
    nvecs = length(synthConfig.link(ll).pt_nums);
    delta = zeros(3,size(Stereo(1).pts,2),nvecs,npair);
    for vec = 1:nvecs
        for pair = 1:npair
            delta(:,:,vec,pair) = Stereo(pair).pts(:,:,synthConfig.link(ll).pt_nums(vec))-Stereo(pair).pts(:,:,synthConfig.link(ll).pt_nums(1));
        end
    end
    mean_delta = nanmean(delta,4);
    std_delta_pair = nanstd(delta,0,4);
    std_delta_time = nanstd(delta,0,2);
    mean_delta = nanmean(mean_delta,2);
    mean_delta = squeeze(mean_delta);

    %align BF basis with body points
    if ll ==1 
        z_hat = cross(mean_delta(:,2),mean_delta(:,3))/norm(cross(mean_delta(:,2),mean_delta(:,3)));
        b     = (mean_delta(:,2)+mean_delta(:,3))/2;
        y_hat = cross(z_hat,b)/norm(cross(z_hat,b));
        x_hat = cross(y_hat,z_hat);
        vectors = 1*[x_hat,y_hat,z_hat]'*mean_delta;
    elseif ll==2
        vectors(:,1) = [0,0,0]';
        vectors(:,2) = 1.1*[0,-norm(mean_delta(:,2)),0]';
    elseif ll==3
        vectors(:,1) = [0,0,0]';
        vectors(:,2) = 1.1*[-norm(mean_delta(:,2)),0,0]';
        vectors(:,3) = 1.1*[-norm(mean_delta(:,3)),0,0]';
    elseif ll==4
        vectors(:,1) = [0,0,0]';
        vectors(:,2) = 1.3*[-norm(mean_delta(:,2)),0,0]';
        vectors(:,3) = 1.3*[-norm(mean_delta(:,3)),0,0]';
    elseif ll==5
        vectors(:,1) = [0,0,0]';
        vectors(:,2) = 1.3*[-norm(mean_delta(:,2)),0,0]';
        vectors(:,3) = 1.3*[-norm(mean_delta(:,3)),0,0]';
    elseif ll==6
        vectors(:,1) = [0,0,0]';
        vectors(:,2) = 1.3*[-norm(mean_delta(:,2)),0,0]';
        %vectors(:,3) = 1.15*[-norm(mean_delta(:,3)),0,0]';
    end
    
    synthConfig.link(ll).BFvecs = vectors;
end

% synthConfig.link(1).BFvecs(:,4) = [30,0,-5]';
% synthConfig.link(1).BFvecs(:,5) = [30,0,5]';

%link 2
% synthConfig.link(2).BFvecs(:,1) = [0,0,0]';
% synthConfig.link(2).BFvecs(:,2) = [0,-hum_len,0]';
% 
% %link 3
% synthConfig.link(3).BFvecs(:,1) = [0,0,0]';
% synthConfig.link(3).BFvecs(:,2) = [-rad_len,0,0]';
% 
% %link 4
% synthConfig.link(4).BFvecs(:,1) = [0,0,0]';
% synthConfig.link(4).BFvecs(:,2) = [-met3_len,0,0]';
% 
% %link 5
% synthConfig.link(5).BFvecs(:,1) = [0,0,0]';
% synthConfig.link(5).BFvecs(:,2) = [-met4_len,0,0]';
% 
% %link 6
% synthConfig.link(6).BFvecs(:,1) = [0,0,0]';
% synthConfig.link(6).BFvecs(:,2) = [-met5_len,0,0]';
%% Define DH params
nn = 1;
%----------------------------------Base Body CF Defn---------------------------------------
synthConfig.link(nn).thetas  = [0;0;0;0;0;0];
synthConfig.link(nn).alphas  = [pi/2; pi/2; 0; pi/2; pi/2; 0];
synthConfig.link(nn).disps   = [0;0;0;0;0;0];
synthConfig.link(nn).offsets = [0;0;0;0;0;0];
synthConfig.link(nn).H       = DHTransforms(synthConfig.link(nn).thetas,synthConfig.link(nn).alphas,synthConfig.link(nn).disps,synthConfig.link(nn).offsets);
nn = nn+1;

%----------------------------------Humerus CF Defn---------------------------------------
%hum_len = 60;
synthConfig.link(nn).thetas  = [0;0;0];
synthConfig.link(nn).alphas  = [pi/2; -pi/2; pi/2];
synthConfig.link(nn).disps   = [0;0;norm(synthConfig.link(2).BFvecs(:,2))];
synthConfig.link(nn).offsets = [0;0;0];
synthConfig.link(nn).H = DHTransforms(synthConfig.link(nn).thetas,synthConfig.link(nn).alphas,synthConfig.link(nn).disps,synthConfig.link(nn).offsets);
nn = nn+1;
%----------------------------------Raduis CF Defn---------------------------------------
%rad_len = 70;
gamma = 85/180*pi;
synthConfig.link(nn).thetas  = [0];
synthConfig.link(nn).alphas  = [0];
synthConfig.link(nn).disps   = [0];
synthConfig.link(nn).offsets = norm(synthConfig.link(3).BFvecs(:,3));
synthConfig.link(nn).H       = DHTransforms(synthConfig.link(nn).thetas,synthConfig.link(nn).alphas,synthConfig.link(nn).disps,synthConfig.link(nn).offsets);
nn = nn+1;

%----------------------------------Digit 3 Metacarpal CF Defn---------------------------------------
met3_len   = norm(synthConfig.link(4).BFvecs(:,3));
synthConfig.link(nn).thetas  = [0;0];
synthConfig.link(nn).alphas  = [pi/2;0];
synthConfig.link(nn).disps   = [0;0];
synthConfig.link(nn).offsets = [0;met3_len];
synthConfig.link(nn).H = DHTransforms(synthConfig.link(nn).thetas,synthConfig.link(nn).alphas,synthConfig.link(nn).disps,synthConfig.link(nn).offsets);
nn = nn+1;

%----------------------------------Digit 4 Metacarpal CF Defn---------------------------------------
met4_len   = norm(synthConfig.link(5).BFvecs(:,3));
synthConfig.link(nn).thetas  = [0;0];
synthConfig.link(nn).alphas  = [pi/2;0];
synthConfig.link(nn).disps   = [0;0];
synthConfig.link(nn).offsets = [0;met4_len];
synthConfig.link(nn).H = DHTransforms(synthConfig.link(nn).thetas,synthConfig.link(nn).alphas,synthConfig.link(nn).disps,synthConfig.link(nn).offsets);
nn = nn+1;
%----------------------------------Digit 5 Metacarpal CF Defn---------------------------------------
met5_len   = norm(synthConfig.link(6).BFvecs(:,2));
synthConfig.link(nn).thetas  = [0;0];
synthConfig.link(nn).alphas  = [pi/2;0];
synthConfig.link(nn).disps   = [0;0];
synthConfig.link(nn).offsets = [0;met5_len];
synthConfig.link(nn).H = DHTransforms(synthConfig.link(nn).thetas,synthConfig.link(nn).alphas,synthConfig.link(nn).disps,synthConfig.link(nn).offsets);
nn = nn+1;




