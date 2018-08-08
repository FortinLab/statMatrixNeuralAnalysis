function varargout = AssignLFPsiteIDs_OE(varargin)
% ASSIGNLFPSITEIDS_OE MATLAB code for AssignLFPsiteIDs_OE.fig
%      ASSIGNLFPSITEIDS_OE, by itself, creates a new ASSIGNLFPSITEIDS_OE or raises the existing
%      singleton*.
%
%      H = ASSIGNLFPSITEIDS_OE returns the handle to a new ASSIGNLFPSITEIDS_OE or the handle to
%      the existing singleton*.
%
%      ASSIGNLFPSITEIDS_OE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ASSIGNLFPSITEIDS_OE.M with the given input arguments.
%
%      ASSIGNLFPSITEIDS_OE('Property','Value',...) creates a new ASSIGNLFPSITEIDS_OE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before AssignLFPsiteIDs_OE_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to AssignLFPsiteIDs_OE_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help AssignLFPsiteIDs_OE


% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @AssignLFPsiteIDs_OE_OpeningFcn, ...
                   'gui_OutputFcn',  @AssignLFPsiteIDs_OE_OutputFcn, ...
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

%%

%%%---------------------------------------------------------------------%%%
%%-----------------------------------------------------------------------%%
%-------------------------------INITIALIZE--------------------------------%
%%-----------------------------------------------------------------------%%
%%%---------------------------------------------------------------------%%%

%%

function AssignLFPsiteIDs_OE_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);

function ratName_INPUT_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function sessionID_INPUT_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function fileList_LST_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function regionID_INPUT_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function probeNum_INPUT_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function varargout = AssignLFPsiteIDs_OE_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output.UserData;

%%

%%%---------------------------------------------------------------------%%%
%%-----------------------------------------------------------------------%%
%--------------------------------CALLBACKS--------------------------------%
%%-----------------------------------------------------------------------%%
%%%---------------------------------------------------------------------%%%

%%
%% BUTTONS %%
function loadData_BTN_Callback(hObject, eventdata, handles)
handles = guidata(hObject);
fileDir  = uigetdir;
if fileDir==0
    return
else
    files = dir(fileDir);
    fileNames = {files.name};
    contFileLog = cellfun(@(a)~isempty(a), strfind(fileNames, '100_CH'));
    contFiles = fileNames(contFileLog);
    handles.fileList_LST.String = contFiles';
    handles.output.UserData = [contFiles', cell(length(contFiles),2)];
end
[ratNameStart, ratNameEnd] = regexp(fileDir, 'Age1-([0-9]*)');
handles.ratName_INPUT.String = fileDir(ratNameStart:ratNameEnd);
[sessionIDstart, sessionIDend] = regexp(fileDir, 'Session([0-9]*[a-z]*)');
handles.sessionID_INPUT.String = fileDir(sessionIDstart:sessionIDend);
%%% Save the directory pathing so that it can be used later %%%
guidata(hObject, handles);

function removeChan_BTN_Callback(hObject, eventdata, handles)
handles = guidata(hObject);
curRecTags = handles.output.UserData;
curChanList = handles.fileList_LST.String;
curChanNum = handles.fileList_LST.Value;

curRecTags(curChanNum,:) = [];
curChanList(curChanNum) = [];

handles.output.UserData = curRecTags;
handles.fileList_LST.String = curChanList;
if curChanNum == 1
    handles.fileList_LST.Value = 1;
else
    handles.fileList_LST.Value = curChanNum-1;
end
guidata(hObject,handles);

function assignID_BTN_Callback(hObject, eventdata, handles)
handles = guidata(hObject);
curChanNum = handles.fileList_LST.Value;
curRecTags = handles.output.UserData;
curRecTags{curChanNum,2} = handles.regionID_INPUT.String;
curRecTags{curChanNum,3} = handles.probeNum_INPUT.String;
handles.output.UserData = curRecTags;
guidata(hObject,handles);

function startCompile_BTN_Callback(hObject, eventdata, handles)
handles = guidata(hObject);
curRecTags = handles.output.UserData;

% Check for unidentified files
unIDdFilesLog = cellfun(@(a)isempty(a), curRecTags(:,2));
unIDdFiles = curRecTags(unIDdFilesLog,1);

if ~isempty(unIDdFiles)
    d = msgbox('Not all selected channels identified. Please identify all files in the list and remove unused ones.');
    uiwait(d);
else
    %%%% Create statMatrix files here %%%%
    regions = unique(curRecTags(:,2));
    % First check to make sure all the channels are identified with unique
    % region/probe IDs
    for r = 1:length(regions)
        curRegionLog = strcmp(curRecTags(:,2), regions{r});
        chans = curRecTags(curRegionLog,3);
        if length(chans)>length(unique(chans))
            d = msgbox(['Multiple channels identified with the same probe number for Region = ' regions{r}]);
            uiwait(d);
        else
        end
    end
    % Once all regions are checked for uniqe Probe IDs now make the
    % behavior matrix that will be associated with each file
    %% First create all the stuff that's common for all matrices
    %%% Timestamp Vector %%%
    %%% Behavior Matrix %%%
    % Make sure to save the behavior matrix separately as well.
    %% Now create the statMatrix file for each region
    % Now create the statMatrix for all the probe regions
    for r = 1:length(regions)
        curRegionLog = strcmp(curRecTags(:,2), regions{r});
        curRegionFiles = curRecTags(curRegionLog,:);
        chans = curRegionFiles(:,3);
%         tempRegionStatMatrix = zeros(length(tsVect),numberOfColumns);
        for c = 1:length(chans)
            curChanFile = curRegionFiles{c,1};
            curChanID = [curRegionFiles{c,2} curRegionFiles{c,3}];
            % Put in statMatrix creation here
        end
    end
end





%% INPUT FIELDS %%
function ratName_INPUT_Callback(hObject, eventdata, handles)

function sessionID_INPUT_Callback(hObject, eventdata, handles)

function regionID_INPUT_Callback(hObject, eventdata, handles)

function probeNum_INPUT_Callback(hObject, eventdata, handles)

%% LIST %%
function fileList_LST_Callback(hObject, eventdata, handles)
handles = guidata(hObject);
curChanNum = handles.fileList_LST.Value;
curRecTags = handles.output.UserData;
handles.regionID_INPUT.String = curRecTags{curChanNum,2};
handles.probeNum_INPUT.String = curRecTags{curChanNum,3};
guidata(hObject,handles);


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in listbox3.
function listbox3_Callback(hObject, eventdata, handles)
% hObject    handle to listbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox3


% --- Executes during object creation, after setting all properties.
function listbox3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
