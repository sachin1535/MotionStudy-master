%slice tests
tests = [1,2,3,4,5,6,10,11];
for tt = tests
    test_num = num2str(tt);
    while length(test_num)<3
        test_num = ['0',test_num];
    end
    cam_folders = dir(['.\Test',test_num,filesep,'Cam*.']);
    for cc = 1:length(cam_folders)
        fprintf('Delacing Camera %d from Test %d...\n',cc,tt)
        cam_num = cam_folders(cc).name(end-2:end);
        vid_slice(['.',filesep,'Test',test_num,filesep,cam_folders(cc).name],['cam',cam_num,'.MP4'],'png')
    end
end
