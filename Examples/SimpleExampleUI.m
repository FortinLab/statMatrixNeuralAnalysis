function SimpleExampleUI
%%
% uiFig = figure('Name', 'Simple UI', 'toolbar', 'none');
uiFig = figure('Name', 'Simple UI', 'DeleteFcn', @WindowClose);

annotation(uiFig, 'textbox', 'Position',[0.25 0.94 0.5 0.05], 'String', 'Simple UI',...
    'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom',...
    'FitHeightToText','off','LineStyle','none');

uiData.leftButton = uicontrol(uiFig, 'Units', 'Normalized', 'Position', [0.05 0.2 0.4 0.75],...
    'Style', 'pushbutton', 'String', 'Left', 'FontSize',12, 'FontWeight', 'bold',...
    'Callback', @ActivateLeft);
uiData.rightButton = uicontrol(uiFig, 'Units', 'Normalized', 'Position', [0.55 0.2 0.4 0.75],...
    'Style', 'pushbutton', 'String', 'Right', 'FontSize',12, 'FontWeight', 'bold',...
    'Callback', @ActivateRight);

annotation(uiFig, 'textbox', 'position', [0.05 0.14 0.4 0.05], 'String', 'Left Duration',...
    'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom',...
    'FitHeightToText', 'off', 'LineStyle', 'none');
annotation(uiFig, 'textbox', 'position', [0.55 0.14 0.4 0.05], 'String', 'Right Duration',...
    'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom',...
    'FitHeightToText', 'off', 'LineStyle', 'none');
uiData.leftDuration = uicontrol(uiFig, 'Units', 'Normalized', 'Position', [0.05 0.05 0.4 0.1],...
    'Style', 'edit', 'String', '100', 'Callback', @UpdateLeftDuration);
uiData.rightDuration = uicontrol(uiFig, 'Units', 'Normalized', 'Position', [0.55 0.05 0.4 0.1],...
    'Style', 'edit', 'String', '100', 'Callback', @UpdateRightDuration);

% Read up on the timer object in the matlab documentation
uiData.leftRewardTimer = timer('StartFcn', {@WaterCallback, 0, 1}, 'TimerFcn', {@WaterCallback, 0, 0}, 'StartDelay', 0.1);
uiData.rightRewardTimer = timer('StartFcn', {@WaterCallback, 1, 1}, 'TimerFcn', {@WaterCallback, 1, 0}, 'StartDelay', 0.1);

uiData.sessionStart = tic;
uiData.leftRewardTimes = nan(1,10);
uiData.rightRewardTimes = nan(1,10);
uiFig.UserData = uiData;

end

function ActivateLeft(hObject, evt)
fig = get(hObject, 'parent');
handles = get(fig, 'UserData');
start(handles.leftRewardTimer);
lastLR = find(~isnan(handles.leftRewardTimes),1,'last');
if lastLR == size(handles.leftRewardTimes,2)
    handles.leftRewardTimes = [handles.leftRewardTimes, nan(1,10)];
end
handles.leftRewardTimes(lastLR+1) = toc(handles.sessionStart);
end

function ActivateRight(hObject,evt)
fig = get(hObject, 'parent');
handles = get(fig, 'UserData');
start(handles.rightRewardTimer);
end

function UpdateLeftDuration(hObject, evt)
fig = get(hObject, 'parent');
handles = get(fig, 'UserData');
handles.leftRewardTimer.StartDelay = str2double(handles.leftDuration.String)/1000;
end

function UpdateRightDuration(hObject, evt)
fig = get(hObject, 'parent');
handles = get(fig, 'UserData');
handles.rightRewardTimer.StartDelay = str2double(handles.rightDuration.String)/1000;
end

function WaterCallback(obj, event, pinID, setting)
if pinID == 0
    if setting == 1
        % Replace this with code to set Arduino pin for left water high
        disp('Left Water On');
    elseif setting == 0
        % Replace this with code to set Arduino pin for left water low
        disp('Left Water Off');
    end
elseif pinID == 1
    if setting == 1
        % Replace this with code to set Arduino pin for right water high
        disp('Right Water On');
    elseif setting == 0
        % Replace this with code to set Arduino pin for right water low
        disp('Right Water Off');
    end
end
end

function WindowClose(hObject, evt)
% Write code here to save whatever you want from the file.
fig = get(hObject, 'parent');
handles = get(fig, 'UserData');
disp('Window Closed');
end