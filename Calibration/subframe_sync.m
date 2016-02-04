function camstruct = subframe_sync(camstruct)
fs = 120;
%determine the number of cameras in the structure
ncam = length(camstruct);
cams = [];
%make a list of cams which have sync information available
for cc = 1:ncam
    if ~isempty(camstruct(cc).sync_del)
        cams = [cams,cc];
    end
end
%for all cams with sync info available,
for cc = cams
    [~,nsteps, npts] = size(camstruct(cc).pts_rect);        %determine the number of points and number of steps in this camera
    del = camstruct(cc).sync_del*fs;                    %determine how many frames of delay exist
    sub_del = del-floor(del);                               %determine the subframe delay 
    %assbemble a matrix to perform the following computation
    sync_mat = (1-sub_del)*[eye(nsteps-1);zeros(1,nsteps-1)]+[zeros(1,nsteps-1);sub_del*eye(nsteps-1)];
    for pp = 1:npts
        A = camstruct(cc).pts_rect(:,:,pp);
        nanids = isnan(camstruct(cc).pts_rect(:,:,pp));
        sum_nans = sum(nanids);
        num2naninds = sum_nans*([eye(nsteps-1);zeros(1,nsteps-1)]-[zeros(1,nsteps-1);eye(nsteps-1)])<0;
        num2naninds = [num2naninds,0;num2naninds,0];
        A(nanids) = 0; 
        B = [A*sync_mat,[NaN;NaN]];
        B(nanids|num2naninds) = NaN;
        camstruct(cc).pts_sync(:,:,pp) = B(:,1:camstruct(cc).end_frame-camstruct(cc).start_frame+1);
    end
end

