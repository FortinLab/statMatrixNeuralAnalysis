function manualArtifactRemoval
global LFP_trace LFPplot
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

% Plot LFP trace

% figure('units','normalized','outerposition',[0 0 1 1]) %For full a fullscreen figure
figure(1)
LFP_trace = statMatrix(:,2);
LFPplot = plot(LFP_trace);

% UI buttons for manual artifact detection and removal
saveLimitsbtn = uicontrol('Style', 'pushbutton', 'String', 'Get Axes Limits',...
    'Position', [100,360,100,20],'Callback', @saveAxisLimits);
clearLimitsbtn = uicontrol('Style', 'pushbutton', 'String', 'Clear Stored Axes',...
    'Position', [100,330,120,20],'Callback', @clearAxisLimits);
removeBadInxsbtn = uicontrol('Style', 'pushbutton', 'String', 'Remove Unwanted Indices',...
    'Position', [100,300,150,20],'Callback', @removeUnwantedIndx);
updatePlotbtn = uicontrol('Style', 'pushbutton', 'String', 'Refresh Plot',...
    'Position', [100,270,100,20],'Callback', @updatePlot);

% UI callback functions
    function saveAxisLimits(source,event)
        global BadIndx
        currentXLimits = {get(gca,'xlim')};
        BadIndx = vertcat(BadIndx,currentXLimits);
        disp('X-axis limits storage updated')
        disp(currentXLimits{1})
        fprintf('Axes count: %d\n',length(BadIndx))
    end

    function clearAxisLimits(source,event)
        global BadIndx
        BadIndx = [];
        disp('Stored axes limits have been cleared')
    end

    function removeUnwantedIndx(source,event)
        global BadIndx
        BadIndxs = zeros(size(LFP_trace));
        for ndx = 1:length(BadIndx)
            BadIndxVector = BadIndx{ndx};
            % Compensates for the case that the captured x-limits exceed
            % length of the LFP trace
            if BadIndxVector(2)> length(LFP_trace)
                Index_Init = round(BadIndxVector(1));
                Index_Fin = length(LFP_trace);
            else
                Index_Init = round(BadIndxVector(1));
                Index_Fin = round(BadIndxVector(2));
            end
            BadIndxs(Index_Init:Index_Fin) = true;
        end
        LFP_trace(BadIndxs==1) = [];
        disp('LFP trace has been updated with removed indices')
    end

    function updatePlot(source,event)
        set(LFPplot, 'YData', LFP_trace);
    end

end

