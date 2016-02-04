function meas = create_meas_matrix(camstruct, options)
link = options.link;
%determine the number of cameras
ncam = length(camstruct);
%determine the number of links to be estimated
nlinks = length(link);
%determine the number of timesteps 
nsteps =  options.tstop - options.tstart + options.interp;
%seed a measurement matrix
meas = zeros(ncam*options.nmeas,nsteps);

for cc = 1:ncam                                             %loop over the cameras
    for ll = 1:nlinks                                       %loop over the links
        npts = size(link(ll).BFvecs,2);                        %determine the number of points on this link
        for pp = 1:npts                                     %loop over the points on this link
        link_inds = camstruct(cc).pt_assoc(1,:) == ll; %find the link associations in this camera
        point_inds = camstruct(cc).pt_assoc(2,:) == pp;%find the point associations in this camera
        [~, point] = find(point_inds & link_inds);          %find where the link and points associations intersect
        ii = 2*(pp-1)+1;                                    %MeasInds to grab
            if isempty(point)                               %This point on this link was not seen, fill with NaN
                keyboard
                indx = options.nmeas*(cc-1)+link(ll).MeasInds(ii:ii+1);
                meas(indx,:) = NaN*ones(2,nsteps);
            else                                            %Otherwize, copy the points
                indx = options.nmeas*(cc-1)+link(ll).MeasInds(ii:ii+1);
                meas(indx,:) = camstruct(cc).pts(:,:,point); 
            end
        end
    end
end


