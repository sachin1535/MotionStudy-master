function Cam = BatIsolation_PointRecognition_local_mb(Cam, options)
%% BatIsolation_PointRcognition
% BatIsolation_PointRcognition.m does the pre-processing and the feature
% extraction of frame images of bats flight motion. The script consists of
% three different parts: the initialization, the image pre-processing and
% the feature extraction.
% Funtion: Pre-processing and Feature Extraction
% *Note that* the configuration in image processing is based on the specific
% condition of bat flight motion image data. It is not assumed to be fitted
% to any other image data. User who wants to use this script should
% consider re-config the parameters of pre-processing part and even some
% parameters in Feature Extraction.

%% 	Initialization
% Initialization part does the basic configuration of the script, the path,
% the folder to read for example. User who use the same archive as the
% author's do not have to change this part but it would be necessary to
% modify if the archive is not same.

working = options.working;
datadir = options.path;

CamNo = options.cams;
  % Last frame you want to process

for cc = CamNo
    SFrameNo=Cam(cc).start_frame; % First frame you want to process
    EFrameNo=Cam(cc).end_frame;
    cam_str = num2str(cc);
    while length(cam_str)<3
        cam_str = ['0',cam_str];
    end
    
    bg=imread([datadir,filesep,'Cam',cam_str,filesep,'bckgnd.png']);% Back ground image
    bg=bg(:,:,1);
    h1=fspecial('gaussian',100,30); % Gaussian filter configuration
    bg1=imfilter(bg,h1); % Back ground image blurred
    
    for FrameNo =  SFrameNo : EFrameNo - 1
        %% Image Pre-processing
        % The image pre-processing part processes the images from original
        % RGB*3 images into *black and white* images contains only the Feature
        % points.
        %cd (['Cam',num2str(cc)])
        strnum = num2str(FrameNo);
        if length(strnum)<3;
            strnum = ['0',strnum];
        end
        
        bat=imread([datadir,filesep,'Cam',cam_str,filesep,strnum,'.png']); % Bat image
        bat=bat(:,:,1);
        %figure;imshow(bat);

        %cd ../

        h2=fspecial('gaussian',100,30);
        bat1=imfilter(bat,h2); % Bat image blurred

        bat_nbg = bg1-bat1;

        threshold = graythresh(bat_nbg);
        bwbat0 = im2bw(bat_nbg,threshold); % Black and white pattern of bat's body



        %%
        % Local adaptive threshold segmentation
        %cd('bwbat')
        if exist([datadir,filesep,'Cam',cam_str,filesep,'bwbat',strnum,'.mat'], 'file')
            bwbat = importdata([datadir,filesep,'Cam',cam_str,filesep,'bwbat',strnum,'.mat']);
            RightMax=FarRight(bwbat); % Find out the rightmost pixel of bat's body
            LeftMin=FarLeft(bwbat); % Find out the leftmost pixel of bat's body

        else
            RightMax = Cam(cc).b_box(FrameNo-SFrameNo+1,1)+Cam(cc).b_box(FrameNo-SFrameNo+1,3);
            LeftMin  = Cam(cc).b_box(FrameNo-SFrameNo+1,1);
            bwbat    = zeros(720,1280);
            bwbat(Cam(cc).b_box(FrameNo-SFrameNo+1,2):Cam(cc).b_box(FrameNo-SFrameNo+1,2)+Cam(cc).b_box(FrameNo-SFrameNo+1,4),...
                Cam(cc).b_box(FrameNo-SFrameNo+1,1):Cam(cc).b_box(FrameNo-SFrameNo+1,1)+Cam(cc).b_box(FrameNo-SFrameNo+1,3)) = 1;
        end
        %cd ../
        %cd ../

%         if RightMax<930&&LeftMin>616  % When bat flies into the lighting range:616<y<930
%             % First part
%             [u,M,s]=MeanGraylevel(bwbat,bat,258,1,RightMax,LeftMin); % Calculate average u and M of first part
%             [uu,MM,TT] = MeanPart(bwbat,bat,258,1,u,M,s); 
%             Tm=min(TT(find(TT~=0))); % Choose the minimum threshold
% 
%             % Second part
%             [u1,M1,s1]=MeanGraylevel(bwbat,bat,358,258,RightMax,LeftMin); % Calculate average u and M of second part
%             [uu1,MM1,TT1]=MeanPart(bwbat,bat,358,258,u1,M1,s1);
%             Tm1=min(TT1(find(TT1~=0))); % Choose the minimum threshold
% 
%             % Third part
%             [u2,M2,s2]=MeanGraylevel(bwbat,bat,408,358,RightMax,LeftMin); % Calculate average u and M of third part
%             [uu2,MM2,TT2]=MeanPart(bwbat,bat,408,358,u2,M2,s2);
%             Tm2=min(TT2(find(TT2~=0))); % Choose the minimum threshold
% 
%             % Fourth part
%             [u3,M3,s3]=MeanGraylevel(bwbat,bat,720,408,RightMax,LeftMin); % Calculate u and M of fourth part
%             [uu3,MM3,TT3]=MeanPart(bwbat,bat,713,408,u3,M3,s3);
%             Tm3=min(TT3(find(TT3~=0))); % Choose the minimum threshold
% 
%             % thresholding
%             newbat=zeros(720,1280);
%             if u<8.5&&s<8.5
%                 [newbat]=IM2BW(newbat,bwbat,bat,258,1,RightMax,LeftMin,Tm+2*u);
%             else
%                 [newbat]=IM2BW(newbat,bwbat,bat,258,1,RightMax,LeftMin,(Tm+(u+s))/2);
%             end
%             if u1<8.5&&s1<8.5
%                 [newbat]=IM2BW(newbat,bwbat,bat,358,258,RightMax,LeftMin,Tm1+2*u1);
%             else
%                 [newbat]=IM2BW(newbat,bwbat,bat,358,258,RightMax,LeftMin,(Tm1+(u1+s1)/2)/2);
%             end
%             if u2<8.5&&s2<8.5
%                 [newbat]=IM2BW(newbat,bwbat,bat,408,358,RightMax,LeftMin,Tm2+2*u2);
%             else
%                 [newbat]=IM2BW(newbat,bwbat,bat,408,358,RightMax,LeftMin,(Tm2+(u2+s2)/2)/2);
%             end
%             if u3<8.5&&s3<8.5
%                 [newbat]=IM2BW(newbat,bwbat,bat,713,408,RightMax,LeftMin,Tm3+2*u3);
%             else
%                 [newbat]=IM2BW(newbat,bwbat,bat,713,408,RightMax,LeftMin,(Tm3+(u3+s3))/2);
%             end
% 
%             % de-noise
%             bwbat=double(bwbat);
% 
%             BW=bwareaopen_wjz(newbat,50,26); % Cross off the noises of white blocks whose area is larger than 100 (8 neighbor)
%             BW1=bwareaopen(BW,3,4);
% 
%             %cd ('SeniorDesign_updated_5_22(in use)')
%         else
            %cd ('SeniorDesign_updated_5_22(in use)')
            %%
%             for i = 1 :720
%                 for k = 1 : 1280
%                     bb(i,k)=double(bwbat0(i,k))*bat(i,k);
%                 end
%             end 
            bb=bwbat0.*double(bat);
            [em,level]=otsu_test1(bg,bat);
            %a=0.1*em+(level-0.3647)/20;
            a=0.1*em+(level-0.1)/20;
            newbat=im2bw(bb,a);  % Bat's body with noises

            % de-noise
            bwbat=double(bwbat);
            BW=bwareaopen_wjz(newbat,50,26);
            BW1=bwareaopen(BW,2,4);
        %end

%         for i=1:720
%             F=find(bwbat(i,:)~=0, 1, 'first');
%             L=find(bwbat(i,:)~=0, 1, 'last' );
%             if ~(isempty(F)&&isempty(L))
%                 for j=1:1280
%                     if (j<F+3&&BW1(i,j)==1)||(BW1(i,j)==1&&j>L-4)
%                         BW1(i,j)=0;
%                     end
%                 end
%             end
%         end
        BW1(~bwbat)=0;


        %% Feature Extraction
        % The Feature Extraction is implemented by computing the centroid of
        % white blocks which stands for feature points.

        Points=regionprops(BW1,'Centroid'); % Coordinates of feature points extracted
        l=length(Points);
        figure(1)
        imshow(bat);
        hold on
        % Show the feature points on bat's body
        centroids = reshape([Points.Centroid],2,[]);
        plot(centroids(1,:)',centroids(2,:)','*r');hold on
        % put the feature points in the camera structure
        Cam(cc).im_feat{FrameNo} = centroids;

        %     cd ../
    %     cd ('YS_PointData\LOCAL')
    %     fid=fopen(['Points',num2str(FrameNo)],'w');
    %     
    %     for h = 1 : l % Save the coordinates into txt files
    %         temppoint=Points(h).Centroid;
    %         xx=temppoint(1,1);
    %         yy=temppoint(1,2);
    %         lineh=horzcat('Markers_',num2str(h),' , ',num2str(xx),' , ',num2str(yy));
    %         
    %         fprintf(fid,'%s\r\n' ,lineh);
    %     end
    %     
    %     fclose(fid);
    %     cd ../
    %     cd ../
    %     cd ../
    end
    
end
fprintf('Done Identifying Features.\n')
end
function [RightMax]=FarRight(bwbat)
k=1;
for i = 1:720
    for j = 1:1279
        if bwbat(i,j)==1&&bwbat(i,j+1)==0
            Right(k)=j;
        end
    end
    k=k+1;
end
RightMax=max(Right);
end

function [LeftMin]=FarLeft(bwbat)
k=1;
for i = 1:720
    for j = 2:1280
        if bwbat(i,j)==1&&bwbat(i,j-1)==0
            Left(k)=j;
        end
    end
    k=k+1;
end
Left(Left==0)=NaN;
LeftMin=min(Left);
end

% FrameNo=input('Please input the first frame number you want to process:'); 
% cd ('FrameFolder')

% bg=imread('Frames_1.jpg');% Back ground image
% bg=bg(:,:,1);
%     imshow(bg);
%     title('Background image')
function [em,level]=otsu_test1(bg,bat)
h1=fspecial('gaussian',10,5); % Gaussian filter configuration
bg1=imfilter(bg,h1); % Back ground image blurred

% bat=imread('Frames_302.jpg'); % Bat image
% bat=bat(:,:,1);
%     figure;imshow(bat);
%     title(['Bat_image_',num2str(FrameNo)]);

% cd ../

h2=fspecial('gaussian',10,5);
bat1=imfilter(bat,h2); % Bat image blurred
%     figure;imshow(bat1);
bat_nbg1 = bg1-bat1;
%     bat_nbg = bg-bat;
%     imshow(bat_nbg);% Bat body without background
%     title('No backgroud');

threshold = graythresh(bat_nbg1);
bwbat1=im2bw(bat_nbg1,threshold); % Black and white pattern of bat's body
% test=imread('untitled.jpg');
% for i = 1 :720
%     for k = 1 : 1280
%         bb(i,k)=double(bwbat(i,k))*bat(i,k);
%     end
% end

[wid,len]=size(bat);     %wid为行数，len为列数
% colorlevel=256;     
num_bins = 256;
hist=zeros(num_bins,1);   %histogram

%threshold=128; %初始阈值
%计算直方图，统计灰度值的个数

for i=1:720  %col
    for j=1:1280  %vol
        if bwbat1(i,j)==1
            m=double(bat(i,j))+1;   % Calculate the gray-level value
            hist(m)=hist(m)+1;
        end
    end
end

p = hist / sum(hist);
omega = cumsum(p);
mu = cumsum(p .* (1:num_bins)');
mu_t = mu(end);

sigma_b_squared = (mu_t * omega - mu).^2 ./ (omega .* (1 - omega));
maxval = max(sigma_b_squared);
isfinite_maxval = isfinite(maxval);
if isfinite_maxval
    idx = mean(find(sigma_b_squared == maxval));
    % Normalize the threshold to the range [0, 1].
    level = (idx - 1) / (num_bins - 1);
else
    level = 0.0;
end

if nargout > 1
    if isfinite_maxval
        em = maxval/(sum(p.*((1:num_bins).^2)') - mu_t^2);
    else
        em = 0;
    end
end
end

function Cam = Get_bwbat(Cam, options)

cams = options.cams;

for cc = cams;

fprintf('Creating BW mask for Cam %d...\n',cc)
SFrameNo=Cam(cc).start_frame; % First frame you want to process
EFrameNo=Cam(cc).end_frame;  % Last frame you want to process

cam_num_str = num2str(cc);
while length(cam_num_str)<3
    cam_num_str = ['0',cam_num_str];
end

bg=imread([options.path,filesep,'Cam',cam_num_str,filesep,'bckgnd.jpg']);% Back ground image
bg=bg(:,:,1);

h1=fspecial('gaussian',12,10); % Gaussian filter configuration
bg1=imfilter(bg,h1); % Back ground image blurred
%%
% specify and positions of all cameras
% turn the black camera area into white for latter compensation

h11=fspecial('gaussian',20,20);
bg2=imfilter(bg,h11);
bwbg = ~im2bw(bg2,0.15);

% for i=1:720
%     for j=1:1280
%         if bwbg(i,j)==0
%             bwbg(i,j)=1;
%         else
%             bwbg(i,j)=0;
%         end
%     end
% end
centers=regionprops(bwbg,'Centroid');
pixels=regionprops(bwbg,'PixelList');
%%

for FrameNo =  SFrameNo : EFrameNo
    frame_str = num2str(FrameNo);
    while length(frame_str)<3
        frame_str = ['0',frame_str];
    end
    
    if exist([options.path,filesep,'Cam',cam_num_str,filesep,'bwbat',frame_str,'.mat'],'file')
        fprintf(['BW mask for Cam',cam_num_str,' frame ',frame_str,' already exists.\n'])
        continue
    end
    %% Image Pre-processing
    % The image pre-processing part processes the images from original
    % RGB*3 images into *black and white* images contains only the Feature
    % points.

    
    bat=imread([options.path,filesep,'Cam',cam_num_str,filesep,frame_str,'.jpg']); % Bat image
    bat=rgb2gray(bat);
    %cd ../
    
    h2=fspecial('gaussian',12,10);
    bat1=imfilter(bat,h2); % Bat image blurred
    
    bat_nbg = bg1-bat1;
    
    
    threshold = graythresh(bat_nbg);
    bwbat=im2bw(bat_nbg,threshold);
    
    center=regionprops(bwbat,'Centroid');
    pixel=regionprops(bwbat,'PixelList');
    %%
    % This part is camera compensation
    if ~isempty(center)
        if length(center)>1
            a=length(pixel(1).PixelList);
            I=1;
            for i=1:length(center)
                a=length(pixel(i).PixelList);
                if a<length(pixel(i).PixelList)
                    a=length(pixel(i).PixelList);
                    I=i;
                end
            end
        else
            I=1;
        end
        j=1;n=1;

        for i=1:length(centers)
            dis(i)=sqrt((center(I).Centroid(1,1)-centers(i).Centroid(1,1)).^2+(center(I).Centroid(1,2)-centers(i).Centroid(1,2)).^2);
            if dis(i)<200&&pixels(i).PixelList(1)>center(I).Centroid(1)-50&&pixels(i).PixelList(length(pixels(i).PixelList))<center(I).Centroid(1)+50
                k(j)=i;
                j=j+1;
            elseif dis(i)>=40&&dis(i)<150
                t(n)=i;
                n=n+1;
            end
        end
        if exist('k')
            for p=1:length(k)
                for q=1:length(pixels(k(p)).PixelList)
                    bwbat(pixels(k(p)).PixelList(q,2),pixels(k(p)).PixelList(q,1))=1;
                end
            end
        end
        l=1;
        if exist('t')
            for p=1:length(t)

                o=1;
                % 
                if centers(t(p)).Centroid(1,1)<center(I).Centroid(1,1)
                    Ma=max(pixels(t(p)).PixelList(:,2));
                    Mi=min(pixels(t(p)).PixelList(:,2));
                    for r=Mi:Ma
                        if ~isempty(find(bwbat(r,:)~=0, 1, 'first'))&&~isempty((find(bwbat(r,:)~=0, 1, 'last')))
                            F(o)=find(bwbat(r,:)~=0, 1, 'first');
                            L(o)=find(bwbat(r,:)~=0, 1, 'last');
                            o=o+1;
                        end
                    end
                    if exist('F')&&exist('L')
                        slope=(Mi-Ma)/(F(1)-F(length(F)));
                        co=Mi-slope*F(1);
                        for r=Mi:Ma
                            for rr= min(F):max(F)
                                if slope*rr+co<r-3
                                    bwbat(r,rr)=1;
                                end
                            end
                        end
                    end
                elseif centers(t(p)).Centroid(1,1)>center(I).Centroid(1,1) 
                    Ma=max(pixels(t(p)).PixelList(:,2));
                    Mi=min(pixels(t(p)).PixelList(:,2));
                    for r=Mi:Ma
                        if ~isempty(find(bwbat(r,:)~=0, 1, 'first'))&&~isempty((find(bwbat(r,:)~=0, 1, 'last')))
                            F(o)=find(bwbat(r,:)~=0, 1, 'first');
                            L(o)=find(bwbat(r,:)~=0, 1, 'last');
                            o=o+1;
                        end
                    end
                    if exist('F')&&exist('L')
                        slope=(Mi-Ma)/(L(1)-L(length(L)));
                        co=Mi-slope*L(1);
                        for r=Mi:Ma
                            for rr= min(L):max(L)
                                if slope*rr+co>r+3
                                    bwbat(r,rr)=1;
                                end
                            end
                        end
                    end
                end
            end
        end
        for r=1:720
            if ~isempty(find(bwbat(r,:)~=0, 1, 'first'))&&~isempty((find(bwbat(r,:)~=0, 1, 'last')))
                ff=find(bwbat(r,:)~=0, 1, 'first');
                ll=find(bwbat(r,:)~=0, 1, 'last');
            end
            if exist('ff')&&exist('ll')
                for rr=ff:ll
                    if bwbat(r,rr)~=1
                        bwbat(r,rr)=1;
                    end
                end
            end
            clear('ff');
            clear('ll');
        end
    else
        bwbat = zeros(720,1280);
    end
 %%       
    %figure; imshow(bwbat);
    Name=['bwbat',num2str(FrameNo)];
    %cd('bwbat')
    save([options.path,filesep,'Cam',cam_num_str,filesep,Name],'bwbat');
    %cd ../
    %Cam(cc).bw_mask{FrameNo} = bwbat;
    clear('k');
    clear('t');
    clear('F');
    clear('L');
end
end
fprintf('BW Masks Created.\n')
end
function [newbat]=IM2BW(newbat,bwbat,bat,Highest,Lowest,RightMax,LeftMin,Tm)
for i=Lowest:Highest
    for j=LeftMin:RightMax
        if bwbat(i,j)==1&&bat(i,j)>Tm
            newbat(i,j)=255;
        end
    end
end
end

function bw2 = bwareaopen_wjz(varargin)
%BWAREAOPEN Remove small objects from binary image.
%   BW2 = BWAREAOPEN(BW,P) removes from a binary image all connected
%   components (objects) that have fewer than P pixels, producing another
%   binary image BW2.  This operation is known as an area opening.  The
%   default connectivity is 8 for two dimensions, 26 for three dimensions,
%   and CONNDEF(NDIMS(BW),'maximal') for higher dimensions. 
%
%   BW2 = BWAREAOPEN(BW,P,CONN) specifies the desired connectivity.  CONN
%   may have the following scalar values:  
%
%       4     two-dimensional four-connected neighborhood
%       8     two-dimensional eight-connected neighborhood
%       6     three-dimensional six-connected neighborhood
%       18    three-dimensional 18-connected neighborhood
%       26    three-dimensional 26-connected neighborhood
%
%   Connectivity may be defined in a more general way for any dimension by
%   using for CONN a 3-by-3-by- ... -by-3 matrix of 0s and 1s.  The
%   1-valued elements define neighborhood locations relative to the center
%   element of CONN.  CONN must be symmetric about its center element.
%
%   Class Support
%   -------------
%   BW can be a logical or numeric array of any dimension, and it must be
%   nonsparse.
%
%   BW2 is logical.
%
%   Example
%   -------
%   Remove all objects in the image text.png containing fewer than 50
%   pixels.
%
%       BW = imread('text.png');
%       BW2 = bwareaopen(BW,50);
%       imshow(BW);
%       figure, imshow(BW2)
%
%   See also BWCONNCOMP, CONNDEF, REGIONPROPS.

%   Copyright 1993-2011 The MathWorks, Inc.
%   $Revision: 1.10.4.9 $  $Date: 2011/11/09 16:48:52 $

% Input/output specs
% ------------------
% BW:    N-D real full matrix
%        any numeric class
%        sparse not allowed
%        anything that's not logical is converted first using
%          bw = BW ~= 0
%        Empty ok
%        Inf's ok, treated as 1
%        NaN's ok, treated as 1
%
% P:     double scalar
%        nonnegative integer
%
% CONN:  connectivity
%
% BW2:   logical, same size as BW
%        contains only 0s and 1s.

[bw,p,conn] = parse_inputs(varargin{:});

CC = bwconncomp(bw,conn);
area = cellfun(@numel, CC.PixelIdxList);

idxToKeep = CC.PixelIdxList(area <= p);
idxToKeep = vertcat(idxToKeep{:});

bw2 = false(size(bw));
bw2(idxToKeep) = true;
end

%%%
%%% parse_inputs
%%%
function [bw,p,conn] = parse_inputs(varargin)

narginchk(2,3)

bw = varargin{1};
validateattributes(bw,{'numeric' 'logical'},{'nonsparse'},mfilename,'BW',1);
if ~islogical(bw)
    bw = bw ~= 0;
end

p = varargin{2};
validateattributes(p,{'double'},{'scalar' 'integer' 'nonnegative'},...
    mfilename,'P',2);

if (nargin >= 3)
    conn = varargin{3};
else
    conn = conndef(ndims(bw),'maximal');
end
iptcheckconn(conn,mfilename,'CONN',3)
end

function [u,M,s]=MeanGraylevel(bwbat,bat,Highest,Lowest,RightMax,LeftMin)
hist=zeros(256,1);
Sum=0;
 
for i=Lowest:Highest  %col
%     LeftMin=find(bwbat(i,:)~=0, 1, 'first');
%     RightMax=find(bwbat(i,:)~=0, 1, 'last' );
    if ~(isempty(LeftMin)&&isempty(RightMax))
        for j=LeftMin:RightMax  %vol
        %if bwbat(i,j)==1&&j+5<RightMax
            m=double(bat(i,j+5))+1;    % Calculate the gray-level value
            hist(m)=hist(m)+1;
        %end
        end
    end
end
u=(1:256) * hist/sum(hist);
M=find(hist~=0, 1, 'last' ); % Find out the maximum value of M
for i=Lowest:Highest  %col
    for j=LeftMin:RightMax  %vol
        if bwbat(i,j)==1
            Sum=(double(bat(i,j))-u)^2+Sum;
        end
    end
end
s=sqrt(Sum/sum(hist));
end

function [uu,MM,TT]=MeanPart(bwbat,bat,Highest,Lowest,u,M,s)
hist=zeros(256,1);
Sum=0;k=1;
for i=Lowest:Highest  %col
    LeftMin=find(bwbat(i,:)~=0, 1, 'first');
    RightMax=find(bwbat(i,:)~=0, 1, 'last' );
    if ~(isempty(LeftMin)&&isempty(RightMax))
    for j=LeftMin:RightMax  %vol
        for p=1:8
            for q=1:8
                if bwbat(i+p-1,j+q-1)~=0&&(j+q+5)<RightMax
                    m=double(bat(i+p-1,j+q+5))+1;    %因为灰度为0-255所以+1
                    hist(m)=hist(m)+1;
                end
            end
        end
        w=sum(hist);
        ut=(1:256) * hist;
        uu=ut/w;
        MM=find(hist~=0, 1, 'last' );
        if uu>u&&MM>M/5
            for p=1:8
                for q=1:8
                    if bwbat(i+p-1,j+q-1)~=0&&(j+q+5)<RightMax
                        Sum=(double(bat(i+p-1,j+q+5))-uu)^2+Sum;
                    end
                end
            end
            ss(k)=sqrt(Sum/w);
            TT(k)=MM-ss(k);
            k=k+1;
        end
    end
    end
end
TT(k)=0;
if isequal(TT,0)            
    TT=M/5-s;
end
end









