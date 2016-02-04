function CFPlot(H,s)
%CFPLOT             Plots coordinate frames based on input of homogeneous
%transforms.  Stack H transforms on the third dimension.  All transforms
%should be entered in the base frame.  All CFs will plot in the base frame.

dims = size(H);
nframes = dims(3);
axis = s*eye(3,3);
points = zeros(6,3);
origins = zeros(nframes,3);
%figure
hold on

for jj = 1:nframes
    k = 0;
    for ii = 1:3
        points(ii+k,:) = H(1:3,4,jj);
        origins(jj,:)       = H(1:3,4,jj);
        k = k+1;
        points(ii+k,:) = (H(1:3,4,jj)+H(1:3,1:3,jj)*axis(:,ii));
    end
    
    plot3(points(1:2,1),points(1:2,2),points(1:2,3), 'b-', 'LineWidth', 2);
    plot3(points(3:4,1),points(3:4,2),points(3:4,3), 'r-', 'LineWidth', 2);
    plot3(points(5:6,1),points(5:6,2),points(5:6,3), 'g-', 'LineWidth', 2);
    handle = text(origins(jj,1), origins(jj,2)+25, origins(jj,3)-25, strcat('Cam ',num2str(jj)));
       set(handle, 'HorizontalAlignment', 'left', 'FontSize', 14);
end

plot3(origins(:,1),origins(:,2),origins(:,3),'k-.');
