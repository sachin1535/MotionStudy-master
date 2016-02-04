function camstruct = find_bat(cam_num, camstruct, options)
%FIND_BAT creates a rectangle bounding box around the bat in the specified
%CAM_NUM for all available frames.  DIRECT is the directory which must be
%specified.  CAMSTRUCT is the cam structure returned by the functinon and
%contains a new feild for the min and max corners of the bounding box.
%CAM_NUM is a three digit number string which specifies the camera number
%to compute a bounding box for.

direct = options.path;
images = dir([direct,filesep,'Cam',cam_num,filesep,'*.png']);
numFrames = length(images);
img = imread([direct,filesep,'Cam',cam_num,filesep,images(1).name]);
[H,W,~]=size(rgb2gray(img));

hf=waitbar(0,['Computing Bounding Box For Cam ',cam_num,' ... ']);

threshold =15.7; % get threshold
rect = zeros(numFrames,4);     
for p=1:numFrames-1;
    %%   get difference between two frames
    waitbar(p/(numFrames-1),hf,['Computing Bounding Box For Cam ',cam_num,' ... ']);
    img2=double(rgb2gray(imread([direct,filesep,'Cam',cam_num,filesep,images(p).name])));           % p
    img3=double(rgb2gray(imread([direct,filesep,'Cam',cam_num,filesep,images(p+1).name])));      % p+1
    bw1 = abs(img3-img2);
%     contour=zeros(H,W,'uint8');
%     counts=0;        
%     for k=1:H
%         for m=1: W
%             if bw1(k,m) >= threshold
%                 contour(k,m)=255;
%                 counts = counts+1;
%             end
%         end
%     end
    contour = bw1>=threshold;
    counts = sum(sum(contour));
    %contour = cast(255*contour,'uint8');
    bw=bwmorph(contour,'clean',2);
    if p>1
        clf
    end
    %% difference big enough
    figure(1)
    imshow(img2),hold on
    if counts > 100
        [cx,cy]=centeriod(bw);       %  calculate center
        rect(p,:) = find_rect(bw,cx,cy,0);    
        if ~isnan(rect(p,:))
            rectangle('Position',rect(p,:),'Curvature',[0,0],'EdgeColor','r','LineWidth',2)
        end
    end
end
delete(hf)
camstruct(str2double(cam_num)).b_box = rect;



function  [rect]=find_rect(bw,cx,cy,edge)
edge = edge +20;
[row,col]=size(bw);
minx=round(cx-edge); maxx=round(cx+edge);
miny=round(cy-edge); maxy=round(cy+edge);
count=0;

for k=1:row
    for m=1:col
        if  bw(k,m)  
            if ((k >= miny) && (k<= maxy)) && ((m>=minx) && (m<=maxx))
                count=count+1;
                id(1,count)=m;   %x
                id(2,count)=k;    %y
            end
        end
    end
end

if count > 10
    % value get from id  indicates there is no zero
    new_minx = min(id(1,:));  new_maxx=max(id(1,:));
    new_miny = min(id(2,:));  new_maxy=max(id(2,:));
    if new_minx==minx || new_maxx == maxx || new_miny == miny || new_maxy == maxy
        rect = find_rect(bw,cx,cy,edge);
    else
        rect = [ new_minx new_miny new_maxx-new_minx new_maxy-new_miny];
        %  set boundary condition 
        if rect(1,1)+rect(1,3) >  col
            rect(1,3) = col-rect(1,1);
        end
        if rect(1,2) + rect(1,4) >  row
            rect(1,4) = row- rect(1,2);
        end
        if rect(1,3)<10
            rect(1,3) = 10;
        end
        if rect(1,4)<10
            rect(1,4) = 10;
        end
    end
elseif edge > 720
    new_minx = NaN;  new_maxx=NaN;
    new_miny = NaN;  new_maxy=NaN;
    rect = [ new_minx new_miny new_maxx-new_minx new_maxy-new_miny];
else 
    rect = find_rect(bw,cx,cy,edge);
end



function [cx,cy]=centeriod(pic)
 [row,col]=size(pic);
 sumx=0;
 sumy=0;
 count =0;
 for k=1:row
     for m=1:col
         if pic(k,m)
             count = count +1;
             sumx=m+sumx;
             sumy=k+sumy;
         end
     end
 end
 cx=sumx/count;
 cy=sumy/count;

              
