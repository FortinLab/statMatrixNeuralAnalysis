function manualArtifactRemoval
global LFP_trace LFP_plot rms_line RMS_plot ...
    ErrorVector Error_plot statMatrix file path
% Load file directory and desired file
[file,path] = uigetfile('*.mat');
if isequal(file,0)
    disp('User selected Cancel');
else
    disp(['User selected ', fullfile(path,file)]);
    disp('Loading plot.....')
end
currFileName = strcat(path,file);
load(currFileName,'statMatrix')
ErrorVector = nan(size(statMatrix(:,2)));

% Plot LFP trace

% figure('units','normalized','outerposition',[0 0 1 1]) %For full a fullscreen figure
LFP_trace = statMatrix(:,2);
figure(1)
LFP_plot = plot(LFP_trace,'b');
hold on
rms_line = (rms(LFP_trace)*ones(length(LFP_trace),1)) + 2*std(LFP_trace);
RMS_plot = plot(rms_line,'c','LineWidth',1);
hold on
Error_plot = plot(ErrorVector,'r');
title(file, 'Interpreter', 'none');

% UI buttons for manual artifact detection and removal
saveLimitsbtn = uicontrol('Style', 'pushbutton', 'String', 'Get Axes Limits',...
    'Position', [100,360,100,20],'Callback', @saveAxisLimits);
clearLimitsbtn = uicontrol('Style', 'pushbutton', 'String', 'Clear Stored Axes',...
    'Position', [100,330,120,20],'Callback', @clearAxisLimits);
removeBadInxsbtn = uicontrol('Style', 'pushbutton', 'String', 'Remove Unwanted Indices',...
    'Position', [100,300,150,20],'Callback', @removeUnwantedIndx);
updateRMSbtn = uicontrol('Style', 'pushbutton', 'String', 'Refresh RMS Line',...
    'Position', [100,270,100,20],'Callback', @updateRMS);
changeCHbtn = uicontrol('Style', 'pushbutton', 'String', 'Change current channel',...
    'Position', [100,240,150,20],'Callback', @changeCH);
% Should I add a save file 

% UI callback functions
    function saveAxisLimits(source,event) % Saves current x-axis limits based on what the user zooms into
        global BadIndx
        currentXLimits = {get(gca,'xlim')};
        BadIndx = vertcat(BadIndx,currentXLimits);
        disp('X-axis limits storage updated')
        disp(currentXLimits{1})
        fprintf('Index pair count: %d\n',length(BadIndx))
    end

    function clearAxisLimits(source,event) % Clear current storage of saved x-axis limits
        global BadIndx BadIndxs
        BadIndx = [];
        BadIndxs = [];
        disp('Stored axes index limits have been cleared')
    end

    function removeUnwantedIndx(source,event) % Removes sections from original LFP signal based on saved limits
        global BadIndx BadIndxs 
        BadIndxs = zeros(size(LFP_trace));
        for ndx = 1:length(BadIndx)
            BadIndxVector = BadIndx{ndx};
            % Compensates for the case that the captured x-limits exceed
            % length of the LFP trace
            if BadIndxVector(2)> length(LFP_trace)
                Index_Init = round(BadIndxVector(1)); % Values are rounded as needed for index values
                Index_Fin = length(LFP_trace);
            else
                Index_Init = round(BadIndxVector(1));
                Index_Fin = round(BadIndxVector(2));
            end
            BadIndxs(Index_Init:Index_Fin) = true;
        end
        LFP_trace = statMatrix(:,2); % This line is repeated for this function to prevent LFP_trace from 
                                     % overwritten as a result of repeated
                                     % clicks
        
        ErrorVector(find(BadIndxs)) = LFP_trace(find(BadIndxs));
        set (Error_plot,'YData',ErrorVector);
        
        LFP_trace(BadIndxs==1) = nan;
        set( LFP_plot, 'YData', LFP_trace );
        
        disp('LFP trace has been updated with removed indices')
        disp('Removed sections have been highlighted in red')
    end

    function updateRMS(source,event)
        rms_line = rms(LFP_trace(~isnan(LFP_trace)))*ones(length(LFP_trace),1)...
            + 2*nanstd(LFP_trace);
        set( RMS_plot, 'YData', rms_line );
    end

    function changeCH(source,event)
        global BadIndxs
        title('Loading plot.....');
        [file,path] = uigetfile('*.mat');
        if isequal(file,0)
            disp('User selected Cancel');
        else
            disp(['User selected ', fullfile(path,file)]);
            disp('Loading plot.....')
        end
        currFileName = strcat(path,file);
        load(currFileName,'statMatrix')
        LFP_trace = statMatrix(:,2);
        
        ErrorVector(find(BadIndxs)) = LFP_trace(find(BadIndxs));
        set (Error_plot,'YData',ErrorVector);
        
        LFP_trace(BadIndxs==1) = nan;
        set( LFP_plot, 'YData', LFP_trace );
        title(file, 'Interpreter', 'none');
    end

end

% Add a function that allows to switch channels 
% May be better to use 
% Have two vectors: original LFP and a NaN vector
% For all unwanted index ranges, swap values

% Look at Initialize plot script from Gabe
% AssignInBase
% For new files, make that the Y. data source