function [options, camstruct] = audiosync(camstruct,options)

if ~isfield(options,'est.cams')
    options.est.cams = input('please specify the cameras to sync:');
end
%Design a highpass filter to chop out low freq noise
ii = 0;
cam_nums = [];
max_ts = [];
max_amps = [];
D = designfilt('highpassiir', 'FilterOrder', 2,'PassbandFrequency', 1e4, 'SampleRate', 44100);
for cc = options.est.cams                       %cycle through files
    ii = ii+1;
    cam_nums = [cam_nums,cc];
    filename = [options.path,filesep,'Cam',num2str(cc),filesep,'cam',num2str(cc),'.MP4'];
    fprintf('Reading Audio From %s ...\n', filename);
    [sig, fs(ii)] = audioread(filename); %read the audio file
    %sig_filt{ii} = filtfilt(D,sig); %filter the audio with the fither object D
    sig_filt{ii} = abs(sig);
    [amp_max, I_max] = max(sig_filt{ii}(:,1));  %Find the max point in the filtered signal
    max_ts = [max_ts, I_max];                   %append max_ts with the new index
    max_amps = [max_amps, amp_max];             %append max_amps with the new index
end

max_ts = max_ts./fs;                            %divide by sampling frequency to get time

[t_min, I] = min(max_ts);                       %find the last camera to start recording

delta_ts = max_ts - t_min;                      %start sampling at that camera

% if save_dat
%     save([folder,'_asd.mat'], 'delta_ts', 'cam_nums');
% end

for cc = 1:length(options.est.cams)
    camstruct(options.est.cams(cc)).sync_del = delta_ts(cc); 
end

%% Plotting function
plotflag = 1;
if plotflag
%Compute subplot dimensions
ncols = 5;  %set the number of plots per row
rows = ceil(ii/ncols); %divide number of cams by nrows to determine number of columns
partial = rem(ii,ncols); %determine if there is a partial remainder

if rows ==1;            %if there is only one row, the number of columns is just the partial 
    ncols = partial;
    if ncols == 0
        ncols =5;
    end
end

dt = 5000;
figure                  %create a figure
for rr = 1:rows
    if rr<rows || partial ==0         %for all rows that are not the last row 
        for cc = 1:ncols %for each column
            ind = 5*(rr-1)+cc; %compute the subplot index
            subplot(rows,ncols,ind) %create the subfigure
            plot([2*44100:5*44100]'/fs(ind),sig_filt{ind}(2*44100:5*44100,1),max_ts(ind),max_amps(ind),'+r')
            %plot([max_ts(ind)*fs-dt:max_ts(ind)*fs+dt]'/fs,sig_filt{ind}(max_ts(ind)*fs-dt:max_ts(ind)*fs+dt,1),max_ts(ind),max_amps(ind),'+r') %plot sig, filtered sig, and max
            title(['cam',num2str(cam_nums(ind))]);
            xlabel('time (s)')
            ylabel('amplitude')
        end
    else                %if this is the last row there may not be ncols to plot, so plot the remainder that was computed before
        for cc = 1:partial
            ind = 5*(rr-1)+cc;
            subplot(rows,ncols,ind)
            plot([2*44100:5*44100]'/fs(ind),sig_filt{ind}(2*44100:5*44100,1),max_ts(ind),max_amps(ind),'+r')
            %plot([max_ts(ind)*fs-dt:max_ts(ind)*fs+dt]'/fs,sig_filt{ind}(max_ts(ind)*fs-dt:max_ts(ind)*fs+dt,1),max_ts(ind),max_amps(ind),'+r')
            title(['cam',num2str(cam_nums(ind))]);
            xlabel('time (s)')
            ylabel('amplitude')
        end
    end
end
end
           