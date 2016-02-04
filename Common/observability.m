function nmeas = observability(camstruct, t, pts)

nmeas = zeros(2,length(t),length(pts));
ncam = length(camstruct);

for cc = 1:ncam
    start_frame = camstruct(cc).start_frame;
    dt = start_frame-1-floor(camstruct(cc).sync_del*119.88);
    nmeas = nmeas + ~isnan(camstruct(cc).pts_sync(:,t-dt,pts));
end