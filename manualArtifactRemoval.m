function manualArtifactRemoval
%% manualArtifactRemoval
%   Tool for visualizing statMatrix data and removing artifact periods by
%   examining the raw LFP trace relative to that channel's RMS values.
%
%% Ideas/To Do:
% 1) Add in listbox to enable quick switching between traces
% 2) Add in listbox to show a list of the periods removed, add
%   functionality to highlight that region and then enable you to add
%   selected period back into the "good trace" or to "redefine" the period
% 3) Add in output... probably to save copies of the statMatrix files that
%   have been "cleaned" by this script.

%% Define Globals
global goodDataTrace badDataTrace timeVals badIndx badIndxs statMatrix smFile...
    artifactRemovalFig figAxes goodDataPlot badDataPlot rmsPlot pokePlotHandles...
    selectFilebtn smFileList 

%#ok<*NASGU>
%% Load smFile directory and select the desired smFile
[smFile,smPath] = uigetfile('*.mat');
if isequal(smFile,0)
    disp('User selected Cancel');
else
    disp(['User selected ', fullfile(smPath,smFile)]);
    disp('Loading plot.....')
end
cd(smPath);
currFileName = strcat(smPath,smFile);
load(currFileName,'statMatrix')

files = dir(smPath);
fileNames = {files.name};
% Load the behavior matrix file for poke events plots
behMatFile = fileNames{cellfun(@(a)~isempty(a), strfind(fileNames, 'BehaviorMatrix'))};
load([smPath behMatFile]);
% Identify list of all statMatrix files
smFileList = fileNames(cellfun(@(a)~isempty(a), strfind(fileNames, '_SM')))'; 

%% Initialize variables
goodDataTrace = statMatrix(:,2);
assignin('base', 'goodDataTrace', goodDataTrace);
badDataTrace = nan(size(goodDataTrace));
assignin('base', 'badDataTrace', badDataTrace);
pokePlotHandles = [];
badIndx = [];
badIndxs = [];
timeVals = statMatrix(:,1)-statMatrix(1,1);
assignin('base', 'timeVals', timeVals);
rmsThresh = rms(goodDataTrace) + 2*std(abs(goodDataTrace));

% Pull out poke events data.
pokeEventsCol = cellfun(@(a)~isempty(a), strfind(behavMatrixColIDs, 'PokeEvents'));
pokeInTSs = behavMatrix(behavMatrix(:,pokeEventsCol)==1,1) - statMatrix(1,1); %#ok<NODEF>
pokeOutTSs = behavMatrix(behavMatrix(:,pokeEventsCol)==-1,1) - statMatrix(1,1);
pokeVals = [pokeInTSs, pokeOutTSs];


%% Create Figure
% artifactRemovalFig = figure('units','normalized','outerposition',[0 0 1 1]); %For full a fullscreen figure
artifactRemovalFig = figure;
figAxes = axes(artifactRemovalFig, 'position', [0.2, 0.1, 0.7, 0.8]);
goodDataPlot = plot(figAxes, timeVals, goodDataTrace, 'b', 'YDataSource', 'goodDataTrace');
axis tight;
hold on;
badDataPlot = plot(figAxes, timeVals, badDataTrace, 'r', 'YDataSource', 'badDataTrace');
rmsPlot = line(get(figAxes, 'xlim'), repmat(rmsThresh, [1,2]), 'color', 'k', 'LineWidth', 1);
title(smFile, 'Interpreter', 'none');
PlotPokeVals(figAxes, pokeVals, rmsThresh);


%% Define UI buttons
% UI buttons for index control
saveLimitsbtn = uicontrol(artifactRemovalFig, 'Units', 'Normalized', 'Style', 'pushbutton', 'String', 'Select Bad Indices',...
    'Position', [0.025,0.2,0.075,0.035],'Callback', @SaveAxisLimits); 
removeBadInxsbtn = uicontrol(artifactRemovalFig, 'Units', 'Normalized', 'Style', 'pushbutton', 'String', 'Remove Bad Indices',...
    'Position', [0.1,0.2,0.075,0.035],'Callback', @RemoveBadIndx);
clearLimitsbtn = uicontrol(artifactRemovalFig, 'Units', 'Normalized', 'Style', 'pushbutton', 'String', 'Clear Stored Indices',...
    'Position', [0.0625,0.15,0.075,0.035],'Callback', @ClearStoredIndexes);
updateRMSbtn = uicontrol(artifactRemovalFig, 'Units', 'Normalized', 'Style', 'pushbutton', 'String', 'Refresh RMS Line',...
    'Position', [0.025,0.1,0.15,0.035],'Callback', @UpdateRMS);
indicesListbtn = uicontrol(artifactRemovalFig, 'Units', 'Normalized', 'Style', 'listbox', 'String', smFileList,...
    'Position', [0.0375,0.25,0.125,0.28],'Callback', @selectFile);

% UI buttons for file control
selectFilebtn = uicontrol(artifactRemovalFig, 'Units', 'Normalized', 'Style', 'listbox', 'String', smFileList,...
    'Position', [0.0375,0.65,0.125,0.28],'Callback', @selectFile);
changeCHbtn = uicontrol(artifactRemovalFig, 'Units', 'Normalized', 'Style', 'pushbutton', 'String', 'Change to Selected Channel',...
    'Position', [0.025,0.6,0.15,0.035],'Callback', @ChangeCH);
end

%% UI callback functions
function SaveAxisLimits(source,event) % Saves current x-axis limits based on what the user zooms into
    global badIndx figAxes
    currentXLimits = {get(figAxes,'xlim')};
    badIndx = vertcat(badIndx,currentXLimits);
    disp('X-axis limits storage updated')
    disp(currentXLimits{1})
    fprintf('Index pair count: %d\n',length(badIndx))
end

function ClearStoredIndexes(source,event) % Clear current storage of saved x-axis limits
    global badIndx badIndxs goodDataTrace badDataTrace
    badIndx = [];
    badIndxs = [];
    goodDataTrace = nansum([goodDataTrace, badDataTrace],2);
    assignin('base', 'goodDataTrace', goodDataTrace);
    badDataTrace = nan(size(goodDataTrace));
    assignin('base', 'badDataTrace', badDataTrace);
    refreshdata
    drawnow;    
    disp('Stored axes index limits have been cleared')
end

function RemoveBadIndx(source,event) % Removes sections from original LFP signal based on saved limits
    global badIndx badIndxs goodDataTrace badDataTrace statMatrix timeVals
    badIndxs = false(size(goodDataTrace));
    for ndx = 1:length(badIndx)
        % Compensates for the case that the captured x-limits exceed
%         % length of the LFP trace
%         if badIndx{ndx}(2) >= timeVals(end)
%             initIndx = find(badIndx{ndx}(1)>timeVals,1,'last'); % Values are rounded as needed for index values
%             finIndx = length(goodDataTrace);
%         else
        initIndx = find(badIndx{ndx}(1)>timeVals,1,'last');
        if isempty(initIndx)
            initIndx = 1;
        end
        finIndx = find(badIndx{ndx}(2)<timeVals,1,'first');
        if isempty(finIndx)
            finIndx = length(timeVals);
        end

        %         end
        badIndxs(initIndx:finIndx) = true;
    end
    goodDataTrace = statMatrix(:,2); % This line is repeated for this function to prevent goodDataTrace from
    % overwritten as a result of repeated
    % clicks

    badDataTrace(badIndxs) = goodDataTrace(badIndxs);
    assignin('base', 'badDataTrace', badDataTrace);
%     set (badDataPlot,'YData',badDataTrace);

    goodDataTrace(badIndxs) = nan;
    assignin('base', 'goodDataTrace', goodDataTrace);
%     set(dataPlot, 'YData', goodDataTrace);
    refreshdata
    drawnow
    disp('LFP trace has been updated with removed indices')
    disp('Removed sections have been highlighted in red')
end

function UpdateRMS(source,event)
    global rmsPlot goodDataTrace
    newRMSthresh = rms(goodDataTrace(~isnan(goodDataTrace))) + 2*nanstd(abs(goodDataTrace));
    set(rmsPlot, 'YData', repmat(newRMSthresh, [1,2]));
    UpdatePokePlotVals(newRMSthresh);
end

function selectFile (source,event)
global smListIndex selectFilebtn
smListIndex = selectFilebtn.Value;
end

function ChangeCH(source,event)
    global badIndxs goodDataTrace badDataTrace statMatrix smFile...
        smFileList smListIndex
    title('Loading plot.....');
    curFile = smFileList{smListIndex};
    smFile = curFile;
    curPath = strcat(pwd,'\');
    if isequal(smFile,0)
        disp('User selected Cancel');
    else
        disp(['User selected ', fullfile(curPath,curFile)]);
        disp('Loading plot.....')
    end
    currFileName = strcat(curPath,curFile);
    load(currFileName,'statMatrix')
    goodDataTrace = statMatrix(:,2);

    badDataTrace(badIndxs) = goodDataTrace(badIndxs);
    assignin('base', 'badDataTrace', badDataTrace);
%     set (badDataPlot,'YData',badDataTrace);

    goodDataTrace(badIndxs) = nan;
%     set( dataPlot, 'YData', goodDataTrace );
    assignin('base', 'goodDataTrace', goodDataTrace);
    title(smFile, 'Interpreter', 'none');
    refreshdata
    drawnow
end

%% Plotting Functions
function PlotPokeVals(figAxes, pokeVals, rmsThresh)
    global pokePlotHandles
    for poke = 1:size(pokeVals,1)
        pokePlotHandles{poke} = plot(figAxes, pokeVals(poke,:), repmat(rmsThresh+(rmsThresh*0.1), [1,2]), 'linewidth', 1, 'marker', '*', 'color', 'k');
    end
end

function UpdatePokePlotVals(rmsThresh)
    global pokePlotHandles
    for poke = 1:length(pokePlotHandles)
        set(pokePlotHandles{poke}, 'YData', repmat(rmsThresh + (rmsThresh*0.5), [1,2]));
    end
end
%% Glenn notes
% Add a function that allows to switch channels 
% May be better to use 
% Have two vectors: original LFP and a NaN vector
% For all unwanted index ranges, swap values

% Look at Initialize plot script from Gabe
% AssignInBase
% For new files, make that the Y. data source