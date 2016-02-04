function varargout = PointPropagatorV1(varargin)
% POINTPROPAGATORV1 MATLAB code for PointPropagatorV1.fig
%      POINTPROPAGATORV1, by itself, creates a new POINTPROPAGATORV1 or raises the existing
%      singleton*.
%
%      H = POINTPROPAGATORV1 returns the handle to a new POINTPROPAGATORV1 or the handle to
%      the existing singleton*.
%
%      POINTPROPAGATORV1('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in POINTPROPAGATORV1.M with the given input arguments.
%
%      POINTPROPAGATORV1('Property','Value',...) creates a new POINTPROPAGATORV1 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PointPropagatorV1_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PointPropagatorV1_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PointPropagatorV1

% Last Modified by GUIDE v2.5 15-Jan-2016 10:29:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PointPropagatorV1_OpeningFcn, ...
                   'gui_OutputFcn',  @PointPropagatorV1_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before PointPropagatorV1 is made visible.
function PointPropagatorV1_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PointPropagatorV1 (see VARARGIN)

MStudyHandles = getappdata(0,'MStudyHands');
set(findobj(handles.figure1),'Units','Normalized')

%initialize Cam structure with empty fields
handles.Cam = MStudyHandles.Cam;
for cc = 1:length(handles.Cam)
    if ~isfield(handles.Cam(cc), 'pts')
        handles.Cam(cc).pts = [];
    end
end

%store working directory
handles.working = pwd;
%import handles from MotionStudy GUI

handles.options = MStudyHandles.options;

%store defualt data directory 
handles.data_dir = MStudyHandles.options.path;
filepath_Callback(handles.filepath, eventdata, handles);
handles.browse_press = 0;

% UIWAIT makes PointPropagatorV1 wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% add a continuous value change listener
if ~isfield(handles,'hListener')
    handles.hListener = ...
    addlistener(handles.time_slider,'ContinuousValueChange',@time_slider_Callback);
end

set(handles.current_image,'NextPlot','replacechildren');
%set button down function for image window
%set(handles.current_image, 'ButtonDownFcn', {@pick_points,handles});

handles.lastSliderVal = get(handles.time_slider,'Value');
% Choose default command line output for PointPropagatorV1
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = PointPropagatorV1_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure



function filepath_Callback(hObject, eventdata, handles)
% hObject    handle to filepath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

 

set(hObject, 'String', handles.data_dir);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function filepath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filepath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in browse_button.
function browse_button_Callback(hObject, eventdata, handles)
% hObject    handle to browse_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function current_cam_Callback(hObject, eventdata, handles)
% hObject    handle to current_cam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
current_cam = str2double(get(hObject,'String'));
start_frame = handles.Cam(current_cam).start_frame;
end_frame = handles.Cam(current_cam).end_frame;

if ~isfield(handles.Cam(current_cam),'sync_del')
    set(handles.sync_del,'String','NA')
elseif isempty(handles.Cam(current_cam).sync_del)
    set(handles.sync_del,'String','NA')
else
    set(handles.sync_del,'String',num2str(handles.Cam(current_cam).sync_del*119.88))
end

% If the frames have not been imported, then import them.
if isempty(handles.Cam(current_cam).frames)
    cam_str = num2str(current_cam);
    while length(cam_str)<3
        cam_str = ['0',cam_str];
    end
    %vidObj = VideoReader([handles.data_dir,filesep,'Cam',cam_str,filesep,'cam',cam_str,'.MP4']);
    %color_frames = read(vidObj, [start_frame, end_frame]);
    %nframes = size(color_frames,4);
    nframes = handles.Cam(current_cam).end_frame - handles.Cam(current_cam).start_frame + 1;
    for ff = 1:nframes
        handles.Cam(current_cam).frames(:,:,ff) = double(rgb2gray(imread([handles.data_dir,filesep,'Cam',cam_str,filesep,num2str(ff+start_frame-1),'.png'])))/255;
    end
    
    [~,~,nframes] = size(handles.Cam(current_cam).frames); 
    handles.Cam(current_cam).nframes = nframes;
    handles.Cam(current_cam).rot_img = 0;
    if isempty(handles.Cam(current_cam).pts)
        handles.Cam(current_cam).pts = NaN*ones(2,nframes);
    end
end

%if the camera selected has a start_frame greater than the current time
%value, reset the start_frame value and current timestep to admissible values
if handles.Cam(current_cam).start_frame > get(handles.time_slider,'Value')
    set(handles.time_slider,'Value',handles.Cam(current_cam).start_frame+1);
    set(handles.current_timestep,'String',num2str(handles.Cam(current_cam).start_frame+1));
    set(handles.start_frame,'String',num2str(handles.Cam(current_cam).start_frame));
    set(handles.time_slider, 'Min', handles.Cam(current_cam).start_frame);
else %otherwise, just set to new value
    set(handles.time_slider, 'Min', handles.Cam(current_cam).start_frame);
    set(handles.start_frame, 'String', num2str(handles.Cam(current_cam).start_frame));
end

%if the camera selected has an end_frame greater than the current time
%value, reset the end_frame value and current timestep to admissible values
if handles.Cam(current_cam).end_frame < get(handles.time_slider,'Value')
    set(handles.time_slider,'Value',handles.Cam(current_cam).end_frame-1);
    set(handles.current_timestep,'String',num2str(handles.Cam(current_cam).end_frame-1));
    set(handles.end_frame,'String',num2str(handles.Cam(current_cam).end_frame));
    set(handles.time_slider, 'Max', handles.Cam(current_cam).end_frame);
else %otherwise, just set to new value
    set(handles.time_slider, 'Max', handles.Cam(current_cam).end_frame);
    set(handles.end_frame, 'String', num2str(handles.Cam(current_cam).end_frame));
end

set(handles.rotate_image,'Value',handles.Cam(current_cam).rot_img);
guidata(hObject, handles);
clear vidObj
current_timestep_Callback(handles.current_timestep,eventdata,handles);

% Hints: get(hObject,'String') returns contents of current_cam as text
%        str2double(get(hObject,'String')) returns contents of current_cam as a double

% --- Executes during object creation, after setting all properties.
function current_cam_CreateFcn(hObject, eventdata, handles)
% hObject    handle to current_cam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function current_timestep_Callback(hObject, eventdata, handles)
% hObject    handle to current_timestep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

tstep = str2double(get(hObject,'String'));

if tstep < get(handles.time_slider,'Min')
    tstep = get(handles.time_slider,'Min');
    set(hObject,'String', num2str(tstep));
elseif tstep > get(handles.time_slider,'Max')
    tstep = get(handles.time_slider,'Max');
    set(hObject,'String', num2str(tstep));
end

set(handles.time_slider,'Value',tstep);
%current_cam = str2double(get(handles.current_cam,'String'));
%h = imshow(handles.Cam(current_cam).frames(:,:,:,tstep-handles.Cam(current_cam).start_frame+1),'Parent',handles.current_image);
%handles.image_handle = h;
handles = show_image(hObject, handles);
set(handles.image_handle, 'ButtonDownFcn', {@current_image_ButtonDownFcn, handles});
set(handles.image_handle, 'Interruptible', 'on');
%set(handles.image_handle, 'HitTest', 'off');
% set(h, 'ButtonDownFcn', {@pick_points,handles});
% set(handles.current_image,'Children',h);
guidata(hObject,handles);

% Hints: get(hObject,'String') returns contents of current_timestep as text
%        str2double(get(hObject,'String')) returns contents of current_timestep as a double

% --- Executes during object creation, after setting all properties.
function current_timestep_CreateFcn(hObject, eventdata, handles)
% hObject    handle to current_timestep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function current_point_label_Callback(hObject, eventdata, handles)
% hObject    handle to current_point_label (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of current_point_label as text
%        str2double(get(hObject,'String')) returns contents of current_point_label as a double

% --- Executes during object creation, after setting all properties.
function current_point_label_CreateFcn(hObject, eventdata, handles)
% hObject    handle to current_point_label (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in save_button.
function save_button_Callback(hObject, eventdata, handles)
% hObject    handle to save_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% for cc = 1:length(handles.Cam)
%     handles.Cam(cc).frames = [];
% end
MStudyHandles.Cam = handles.Cam;
setappdata(0,'MStudyHands', MStudyHandles);
%varargout{1} = handles.Cam;
%save([handles.data_dir,filesep,'CamStruct.mat'], '-struct', 'handles','Cam','-v7.3')

% --- Executes on slider movement.
function time_slider_Callback(hObject, eventdata, handles)
% hObject    handle to time_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 handles = guidata(hObject);
 % get the slider value and convert it to the nearest integer that is less
 % than this value
 newVal = round(get(hObject,'Value'));
 % set the slider value to this integer which will be in the set {1,2,3,...,12,13}
 set(hObject,'Value',newVal);
 % now only do something in response to the slider movement if the 
 % new value is different from the last slider value
 if newVal ~= handles.lastSliderVal
     % it is different, so we have moved up or down from the previous integer
     % save the new value
     handles.lastSliderVal = newVal;
     guidata(hObject,handles);
    % display the current value of the slider
    %disp(['at slider value ' num2str(get(hObject,'Value'))]);
    set(handles.current_timestep,'String',num2str(newVal))
    current_cam = str2double(get(handles.current_cam,'String'));
    handles = show_image(hObject,handles);
    %set(handles.image_handle, 'ButtonDownFcn', {@pick_points,handles});
    if ~isempty(handles.Cam(current_cam).pts);
        plot_points(handles.image_handle,handles)
    else
        set(handles.image_handle, 'ButtonDownFcn', {@current_image_ButtonDownFcn, handles});
        set(handles.image_handle, 'Interruptible', 'on');
    end
    %set(handles.image_handle, 'HitTest', 'off');
 end
 guidata(hObject, handles);
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function time_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to time_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function start_frame_Callback(hObject, eventdata, handles)
% hObject    handle to start_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of start_frame as text
%        str2double(get(hObject,'String')) returns contents of start_frame as a double
% handles.Cam(str2double(get(handles.current_cam,'String'))).start_frame = str2double(get(hObject,'String'));
% if str2double(get(hObject,'String')) > get(handles.time_slider,'Value')
%     set(handles.time_slider,'Value', str2double(get(hObject,'String')))
%     set(handles.current_timestep,'String',get(hObject,'String'))
% end
% set(handles.time_slider,'Min',str2double(get(hObject,'String')));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function start_frame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to start_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function end_frame_Callback(hObject, eventdata, handles)
% hObject    handle to end_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% if str2double(get(hObject,'String')) < get(handles.time_slider,'Value')
%     set(handles.time_slider,'Value', str2double(get(hObject,'String')))
%     set(handles.current_timestep,'String',get(hObject,'String'))
% end
% 
% if str2double(get(hObject,'String')) <= handles.Cam(str2double(get(handles.current_cam,'String'))).nframes
%     handles.Cam(str2double(get(handles.current_cam,'String'))).end_frame = str2double(get(hObject,'String'));
%     set(handles.time_slider,'Max',str2double(get(hObject,'String')));
% else
%     handles.Cam(str2double(get(handles.current_cam,'String'))).end_frame = handles.Cam(str2double(get(handles.current_cam,'String'))).nframes;
%     set(handles.time_slider,'Max',handles.Cam(str2double(get(handles.current_cam,'String'))).nframes);
% end
guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of end_frame as text
%        str2double(get(hObject,'String')) returns contents of end_frame as a double

% --- Executes during object creation, after setting all properties.
function end_frame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to end_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in rotate_image.
function rotate_image_Callback(hObject, eventdata, handles)
% hObject    handle to rotate_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of rotate_image

current_cam = str2double(get(handles.current_cam,'String'));
handles.Cam(current_cam).rot_img = get(hObject,'Value');
[~,~,~,nframes] = size(handles.Cam(current_cam).frames);

for ii = 1:nframes
    handles.Cam(current_cam).frames(:,:,:,ii) = imrotate(handles.Cam(current_cam).frames(:,:,:,ii),180);
    handles.Cam(current_cam).features{ii}.Locations = ([cos(pi),-sin(pi);sin(pi),cos(pi)]*handles.Cam(current_cam).features{ii}.Locations')';
end

current_time = get(handles.time_slider,'Value');
h = imshow(handles.Cam(current_cam).frames(:,:,:,current_time),'Parent',handles.current_image);
handles.image_handle = h;
%set(handles.image_handle, 'ButtonDownFcn', {@pick_points,handles});
set(handles.image_handle, 'HitTest', 'off');
guidata(hObject, handles);

%--------------------------------------------------------------------------
%                         My Processing Functions
%--------------------------------------------------------------------------
% function pick_points(hObject, eventdata, handles)
% % hObject    handle to rotate_image (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% % PICK_POINTS specifies new points 
% point = get(get(hObject, 'Parent'),'Currentpoint');
% cam = str2double(get(handles.current_cam,'String'));
% timestep = str2double(get(handles.current_timestep,'String'));
% point_num = str2double(get(handles.current_point_label,'String'));
% 
% handles.Cam(cam).pts(:,timestep,point_num) = point(1,1:2)';
% set(handles.current_image,'NextPlot', 'Add');
% ax_child = get(handles.current_image,'Children');
% child_types = get(ax_child,'Type');
% cell_index = regexp('line',child_types);
% 
% if ~isempty(cell_index)
%     index = find([cell_index{:}]);
%     delete(ax_child(index));
% end
% 
% plot(handles.current_image, reshape(handles.Cam(cam).pts(1,timestep,:),1,[]), reshape(handles.Cam(cam).pts(2,timestep,:),1,[]),'ro');
% set(handles.current_image,'NextPlot', 'replacechildren');
% guidata(hObject,handles);

% --- Executes on mouse press over axes background.
function current_image_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to rotate_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% PICK_POINTS specifies new points 
point = get(get(hObject, 'Parent'),'Currentpoint');

%determine camera number, time step, and point number being selected
cam = str2double(get(handles.current_cam,'String'));
timestep = str2double(get(handles.current_timestep,'String'))-handles.Cam(cam).start_frame+1;
point_num = str2double(get(handles.current_point_label,'String'));
%determine number of total steps and number of total points
[~, nsteps, npts] = size(handles.Cam(cam).pts);
%if the zoom feature is being used, compute the appropriate offset
if handles.Cam(cam).zoom
        x = handles.Cam(cam).b_box(timestep,1)-1;
        y = handles.Cam(cam).b_box(timestep,2)-1;
else %otherwise, the offset is nothing
    x = 0;
    y = 0;
end

%if a new point number has been specified, then add page to the point
%matrix
if npts<point_num
   handles.Cam(cam).pts(:,:,(npts+1):point_num) = NaN*ones(2,nsteps,point_num-npts);
end

persistent chk %delcare a persistent variable to determine if a double click has been used
if isempty(chk) %if no immediately previous click
    chk = 1;    %then chk variable is high
    pause(1); %Add a delay to distinguish single click from a double click
    if chk == 1 %if there was only one click
        %compute the appropriate location based on the new click
        locs(:,1) = handles.Cam(cam).features{timestep}.Location(:,1)-x;
        locs(:,2) = handles.Cam(cam).features{timestep}.Location(:,2)-y;
        %compute the distance to the click location 
        dist = locs - repmat(point(1,1:2),size(locs,1),1);
        delta = sum(dist.*dist,2);
        [val,I] = min(delta);
        if val<10^2     %if there is a point within the threshold, use that point
            handles.Cam(cam).pts(:,timestep,point_num) = locs(I,:)'+[x;y];
        end
        plot_points(hObject, handles);
    else
        
    end
    chk = []; %empty the click check
     %plot the new point and label
else %if the register was full, the this was a double click
    chk = 2; %there were two clicks
    handles.Cam(cam).pts(:,timestep,point_num) = point(1,1:2)'+[x;y]; %take the click location 
    plot_points(hObject, handles); %plot points
end

function plot_points(hObject, handles)
cam = str2double(get(handles.current_cam,'String'));
npts = size(handles.Cam(cam).pts,3);
timestep = str2double(get(handles.current_timestep,'String')) - handles.Cam(cam).start_frame+1;

if any(any(any(~isnan(handles.Cam(cam).pts))))
    set(handles.current_image,'NextPlot', 'Add');
    ax_child = get(handles.current_image,'Children');
    child_types = get(ax_child,'Type');

    if iscell(child_types)
    for ii = 1:length(child_types)
        if ~isempty(regexp('line',child_types{ii})) || ~isempty(regexp('text',child_types{ii}));
            delete(ax_child(ii));
        end
    end
    end

    if handles.Cam(cam).zoom
        x = handles.Cam(cam).b_box(timestep,1)-1;
        y = handles.Cam(cam).b_box(timestep,2)-1;
        plot(handles.current_image, reshape(handles.Cam(cam).pts(1,timestep,:)-x,1,[]), reshape(handles.Cam(cam).pts(2,timestep,:)-y,1,[]),'ro');
        for pp = 1:npts
            if ~isnan(handles.Cam(cam).pts(:,timestep,pp))
                h = text(handles.Cam(cam).pts(1,timestep,pp)-x,handles.Cam(cam).pts(2,timestep,pp)-y,num2str(pp), 'Color', 'r', 'Parent', handles.current_image, 'HorizontalAlignment', 'Right');
                set(h, 'HitTest', 'off');
            end
        end
    else
        plot(handles.current_image, reshape(handles.Cam(cam).pts(1,timestep,:),1,[]), reshape(handles.Cam(cam).pts(2,timestep,:),1,[]),'ro');
        for pp = 1:npts
            if ~isnan(handles.Cam(cam).pts(:,timestep,pp))
                h = text(handles.Cam(cam).pts(1,timestep,pp),handles.Cam(cam).pts(2,timestep,pp),num2str(pp), 'Color', 'r', 'Parent', handles.current_image, 'HorizontalAlignment', 'Right');
                set(h, 'HitTest', 'off');
            end
        end
    end

    set(handles.current_image,'NextPlot', 'replacechildren');
    set(handles.image_handle, 'ButtonDownFcn', {@current_image_ButtonDownFcn, handles});
    set(handles.image_handle, 'Interruptible', 'on');
    guidata(handles.figure1,handles);
    
end

function handles = show_image(hObject, handles)
%all of this functionality allows for the user to click points on the
%image.  However, these computations require significant time and cause
%scrolling to lag.  Pursue more efficient format for point picking in the
%future.
cam = str2double(get(handles.current_cam,'String'));
tstep = round(get(handles.time_slider,'Value'))-handles.Cam(cam).start_frame+1;

if isfield(handles.Cam(cam),'b_box')
    b_box = handles.Cam(cam).b_box(tstep,:);
    h = imshow(handles.Cam(cam).frames(b_box(2):b_box(2)+b_box(4),b_box(1):b_box(1)+b_box(3),tstep),'Parent',handles.current_image);
    handles.Cam(cam).zoom = 1;
else
    h = imshow(handles.Cam(cam).frames(:,:,tstep),'Parent',handles.current_image);
end
handles.image_handle = h;

if ~isempty(handles.Cam(cam).pts) && any(any(~isnan(handles.Cam(cam).pts(:,tstep,:))))
    plot_points(handles.image_handle,handles)
else
    set(handles.image_handle, 'ButtonDownFcn', {@current_image_ButtonDownFcn, handles});
    set(handles.image_handle, 'Interruptible', 'on');
end

function cam_list = which_cams(Cams)
%%Determines which cameras have been imported. CAMS is a structure of
%%camera data. The required subfield is CAM.STARTFRAME.  This field is
%%automatically populated when a camera is imported.  CAM_LIST is a row
%%vector of the cams which have been successfully imported.

ncams = length(Cams);
cam_list = [];
%Run loop for length of Cam structure
for cc = 1:ncams
    %If the camera has been imported, store the camera number
    if ~isempty(Cams(cc).startframe)
        cam_list = [cam_list, cc];
    end
end


% --- Executes on button press in track_pts.
function track_pts_Callback(hObject, eventdata, handles)
% hObject    handle to track_pts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Cam = track_points_im(handles.Cam,handles.options.path);
guidata(hObject,handles);


% --- Executes on button press in delete_pt.
function delete_pt_Callback(hObject, eventdata, handles)
% hObject    handle to delete_pt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cam = str2double(get(handles.current_cam,'String'));
timestep = get(handles.time_slider,'value') - handles.Cam(cam).start_frame+1;
pt_num = str2double(get(handles.current_point_label,'String'));

handles.Cam(cam).pts(:,timestep,pt_num) = [NaN;NaN];
plot_points(hObject, handles)
guidata(hObject,handles);


% --- Executes on button press in delete_multi.
function delete_multi_Callback(hObject, eventdata, handles)
% hObject    handle to delete_multi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cam = str2double(get(handles.current_cam,'String'));
pt_nums = input('Which pts would you like to delete?:');
timesteps = input('Which timesteps would you like to delete?:')- handles.Cam(cam).start_frame+1;

handles.Cam(cam).pts(:,timesteps,pt_nums) = NaN*ones(2,length(timesteps),length(pt_nums));
plot_points(hObject, handles)
guidata(hObject,handles);



function sync_del_Callback(hObject, eventdata, handles)
% hObject    handle to sync_del (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sync_del as text
%        str2double(get(hObject,'String')) returns contents of sync_del as a double


% --- Executes during object creation, after setting all properties.
function sync_del_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sync_del (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
