function varargout = AssignLFPsiteIDs2_OE(varargin)
% ASSIGNLFPSITEIDS2_OE MATLAB code for AssignLFPsiteIDs2_OE.fig
%      ASSIGNLFPSITEIDS2_OE, by itself, creates a new ASSIGNLFPSITEIDS2_OE or raises the existing
%      singleton*.
%
%      H = ASSIGNLFPSITEIDS2_OE returns the handle to a new ASSIGNLFPSITEIDS2_OE or the handle to
%      the existing singleton*.
%
%      ASSIGNLFPSITEIDS2_OE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ASSIGNLFPSITEIDS2_OE.M with the given input arguments.
%
%      ASSIGNLFPSITEIDS2_OE('Property','Value',...) creates a new ASSIGNLFPSITEIDS2_OE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before AssignLFPsiteIDs2_OE_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to AssignLFPsiteIDs2_OE_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help AssignLFPsiteIDs2_OE

% Last Modified by GUIDE v2.5 07-Aug-2018 14:40:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @AssignLFPsiteIDs2_OE_OpeningFcn, ...
                   'gui_OutputFcn',  @AssignLFPsiteIDs2_OE_OutputFcn, ...
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

function AssignLFPsiteIDs2_OE_OpeningFcn(hObject, eventdata, handles, varargin)
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

function unIDdChans_LST_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function iDdChans_LST_CreateFcn(hObject, eventdata, handles)
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

function varargout = AssignLFPsiteIDs2_OE_OutputFcn(hObject, eventdata, handles) 
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
    cd(fileDir);
    files = dir(fileDir);
    fileNames = {files.name};
    contFileLog = cellfun(@(a)~isempty(a), strfind(fileNames, '100_CH'));
    contFiles = fileNames(contFileLog);
    handles.unIDdChans_LST.String = contFiles';
    handles.unIDdChans_LST.Value = [];
    handles.iDdChans_LST.String = ' ';
    handles.iDdChans_LST.Value = [];
    handles.output.UserData = [{contFiles'}, {[]}, {fileDir}];
    handles.ratName_INPUT.String = 'Rat Name';
    handles.sessionID_INPUT.String = 'Session ID';
    handles.regionID_INPUT.String = 'Region ID';
    handles.probeNum_INPUT.String = 'Probe #';
%     handles.
end
guidata(hObject, handles);

function removeChan_BTN_Callback(hObject, eventdata, handles)
handles = guidata(hObject);
curUnIDdFls = handles.output.UserData{1};
curChanList = handles.unIDdChans_LST.String;
curChanNum = handles.unIDdChans_LST.Value;

curUnIDdFls(curChanNum) = [];

handles.output.UserData{1} = curUnIDdFls;
handles.unIDdChans_LST.String = handles.output.UserData{1};
if sum(curChanNum == 1) == 1
    handles.unIDdChans_LST.Value = [];
else
    handles.unIDdChans_LST.Value = curChanNum(end)-1;
end
guidata(hObject,handles);

function removeAllChans_BTN_Callback(hObject, eventdata, handles)

function assignID_BTN_Callback(hObject, eventdata, handles)
% Pull handles
handles = guidata(hObject);
% Identify the file of interest
curChanNum = handles.unIDdChans_LST.Value;
curRecFl = handles.output.UserData{1}{curChanNum};
% Add it into the ID'd file cell (output.UserData{2})
handles.output.UserData{2} = sortrows([handles.output.UserData{2}; {curRecFl}, {handles.regionID_INPUT.String}, {handles.probeNum_INPUT.String}], [2 3]);
handles.iDdChans_LST.String = handles.output.UserData{2}(:,1);
handles.iDdChans_LST.Value = find(strcmp(handles.iDdChans_LST.String, curRecFl));
% Now, remove it from the unID'd list
handles.unIDdChans_LST.String(curChanNum) = [];
% As well as the userData{1} list
handles.output.UserData{1}(curChanNum) = [];
% Set the cursor to highlight the file in the new column now.
handles.iDdChans_LST.Value = length(handles.iDdChans_LST.String);
handles.unIDdChans_LST.Value = [];
guidata(hObject,handles);

function removeID_BTN_Callback(hObject, eventdata, handles)
handles = guidata(hObject);
flList = handles.output.UserData{1};
curRecTags = handles.output.UserData{2};
selFlSpot = handles.iDdChans_LST.Value;
handles.iDdChans_LST.Value = [];
handles.iDdChans_LST.String(selFlSpot) = [];
handles.output.UserData{2}(selFlSpot,:) = [];
handles.output.UserData{1} = sortrows([flList; curRecTags(selFlSpot,1)]);
handles.unIDdChans_LST.String = handles.output.UserData{1};
handles.unIDdChans_LST.Value = find(strcmp(handles.unIDdChans_LST.String, curRecTags{1}));

function startCompile_BTN_Callback(hObject, eventdata, handles)
handles = guidata(hObject);
% Check for unidentified files
if ~isempty(handles.output.UserData{1})
    d = msgbox('Not all selected channels identified. Please identify all files in the list and remove unused ones.');
    uiwait(d);
else
    ratName = handles.ratName_INPUT.String;
    ssnID = handles.sessionID_INPUT.String;
    curRecTags = handles.output.UserData{2};
    regionTags = cellfun(@(a,b)[a '_' b], curRecTags(:,2), curRecTags(:,3), 'uniformoutput',0);
    chanMapStruct = struct('FileDir', handles.output.UserData{3},...
        'FileIDs', {[curRecTags(:,1), regionTags]}, 'RatName', ratName,...
        'SessionID', ssnID); %#ok<NASGU>
    uisave('chanMapStruct', [ratName '_Session' ssnID '_ChanMap.mat']);
    d = msgbox('Channel ID map saved!');
    uiwait(d);
end

%% LIST %%
function unIDdChans_LST_Callback(hObject, eventdata, handles)
% Pull handles
handles = guidata(hObject);
% Set the other list value to null
handles.iDdChans_LST.Value = [];
handles.regionID_INPUT.String = [];
handles.probeNum_INPUT.String = [];
guidata(hObject,handles);

function iDdChans_LST_Callback(hObject, eventdata, handles)
% Pull handles
handles = guidata(hObject);
% Set the other list value to null
handles.unIDdChans_LST.Value = [];
% Get the selected value% Get the selected value
curChanNum = handles.iDdChans_LST.Value;
curRecTags = handles.output.UserData{2};
handles.regionID_INPUT.String = curRecTags{curChanNum,2};
handles.probeNum_INPUT.String = curRecTags{curChanNum,3};
