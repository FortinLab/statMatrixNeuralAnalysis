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
    selectFilebtn smFileList storedIndices indicesListbtn smPath

%#ok<*NASGU>
%% Load smFile directory and select the desired smFile
[smFile,smPath] = uigetfile('*.mat');
cd(smPath);
UpdateDataTrace;
files = dir(smPath);
fileNames = {files.name};
% Load the behavior matrix file for poke events plots
behMatFile = fileNames{cellfun(@(a)~isempty(a), strfind(fileNames, 'BehaviorMatrix'))};
load([smPath behMatFile]);
% Identify list of all statMatrix files
smFileList = fileNames(cellfun(@(a)~isempty(a), regexp(fileNames, '_SM\>')))';

%% Initialize variables
pokePlotHandles = [];
badIndx = [];
badIndxs = [];
storedIndices = [];
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
    'Position', [0.025,0.25,0.075,0.035],'Callback', @SaveAxisLimits);
removeBadInxsbtn = uicontrol(artifactRemovalFig, 'Units', 'Normalized', 'Style', 'pushbutton', 'String', 'Remove Bad Indices',...
    'Position', [0.1,0.25,0.075,0.035],'Callback', @RemoveBadIndx);
clearLimitsbtn = uicontrol(artifactRemovalFig, 'Units', 'Normalized', 'Style', 'pushbutton', 'String', 'Clear Stored Indices',...
    'Position', [0.0625,0.15,0.075,0.035],'Callback', @ClearStoredIndexes);
updateRMSbtn = uicontrol(artifactRemovalFig, 'Units', 'Normalized', 'Style', 'pushbutton', 'String', 'Refresh RMS Line',...
    'Position', [0.025,0.1,0.15,0.035],'Callback', @UpdateRMS);
indicesListbtn = uicontrol(artifactRemovalFig, 'Units', 'Normalized', 'Style', 'listbox', 'String', storedIndices,...
    'Position', [0.0375,0.3,0.125,0.28]);
returnIndicesbtn = uicontrol(artifactRemovalFig, 'Units', 'Normalized', 'Style', 'pushbutton', 'String', 'Return Selected Indices',...
    'Position', [0.025,0.2,0.15,0.035],'Callback', @returnIndices);

% UI buttons for file control
selectFilebtn = uicontrol(artifactRemovalFig, 'Units', 'Normalized', 'Style', 'listbox', 'String', smFileList,...
    'Position', [0.0375,0.65,0.125,0.28],'Callback', @selectFile);
changeCHbtn = uicontrol(artifactRemovalFig, 'Units', 'Normalized', 'Style', 'pushbutton', 'String', 'Change to Selected Channel',...
    'Position', [0.025,0.6,0.15,0.035],'Callback', @ChangeCH);
Save2Filebtn = uicontrol(artifactRemovalFig, 'Units', 'Normalized', 'Style', 'pushbutton', 'String', 'Save Current Signal as File',...
    'Position', [0.025,0.05,0.15,0.035],'Callback', @SaveFile);
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
global badIndx badIndxs goodDataTrace badDataTrace storedIndices indicesListbtn
badIndx = [];
badIndxs = [];
storedIndices = [];
set(indicesListbtn,'String',storedIndices)
goodDataTrace = nansum([goodDataTrace, badDataTrace],2);
assignin('base', 'goodDataTrace', goodDataTrace);
badDataTrace = nan(size(goodDataTrace));
assignin('base', 'badDataTrace', badDataTrace);
refreshdata
drawnow;
disp('Stored axes index limits have been cleared')
end

function RemoveBadIndx(source,event) % Removes sections from original LFP signal based on saved limits
global badIndx badIndxs goodDataTrace badDataTrace statMatrix timeVals storedIndices indicesListbtn
for indxPairCount = 1:length(badIndx)
    if indxPairCount == 1
        indxStrings = cellstr(num2str(badIndx{indxPairCount}));
    else
        indxStrings = vertcat(indxStrings,cellstr(num2str(badIndx{indxPairCount})));
    end
end
storedIndices = indxStrings;
set(indicesListbtn,'String',storedIndices)

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

function returnIndices(source,event)
global indicesListbtn badIndx badIndxs goodDataTrace badDataTrace statMatrix timeVals...
    storedIndices
selectedIndxPair = indicesListbtn.Value;
if length(badIndx) >= selectedIndxPair
    initIndx = find(badIndx{selectedIndxPair}(1)>timeVals,1,'last');
    if isempty(initIndx)
        initIndx = 1;
    end
    finIndx = find(badIndx{selectedIndxPair}(2)<timeVals,1,'first');
    if isempty(finIndx)
        finIndx = length(timeVals);
    end
    badIndxs(initIndx:finIndx) = false;
    badIndx(selectedIndxPair) = [];
    storedIndices(selectedIndxPair) = [];
    set(indicesListbtn,'Value',1);
    set(indicesListbtn,'String',storedIndices)
end
origDataTrace = statMatrix(:,2);

badDataTrace(initIndx:finIndx) = nan;
assignin('base', 'badDataTrace', badDataTrace);

goodDataTrace(initIndx:finIndx) = origDataTrace(initIndx:finIndx);
assignin('base', 'goodDataTrace', goodDataTrace);

refreshdata
drawnow
end

function UpdateRMS(source,event)
global rmsPlot goodDataTrace
newRMSthresh = rms(goodDataTrace(~isnan(goodDataTrace))) + 2*nanstd(abs(goodDataTrace));
set(rmsPlot, 'YData', repmat(newRMSthresh, [1,2]));
UpdatePokePlotVals(newRMSthresh);
end

function selectFile(source,event)
global smListIndex selectFilebtn
smListIndex = selectFilebtn.Value;
end

function ChangeCH(source,event)
global badIndxs goodDataTrace badDataTrace smFile...
    smFileList smListIndex figAxes
title(figAxes,'Loading plot.....');
smFile = smFileList{smListIndex};
UpdateDataTrace;
badDataTrace(badIndxs) = goodDataTrace(badIndxs);
assignin('base', 'badDataTrace', badDataTrace);
%     set (badDataPlot,'YData',badDataTrace);

goodDataTrace(badIndxs) = nan;
%     set( dataPlot, 'YData', goodDataTrace );
assignin('base', 'goodDataTrace', goodDataTrace);
title(figAxes, smFile, 'Interpreter', 'none');
refreshdata
drawnow
end

function SaveFile(source,event)
global statMatrix smFile badIndxs

statMatrix_edited = statMatrix;
statMatrix_edited(badIndxs,:) = nan;
outputfilename = strcat('Edited_',smFile);
uisave('statMatrix_edited',outputfilename);


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

function UpdateDataTrace
global smFile goodDataTrace badDataTrace statMatrix
if isequal(smFile,0)
    disp('User selected Cancel');
else
    disp(['User selected ', smFile]);
    disp('Loading plot.....')
end
load(smFile,'statMatrix')

goodDataTrace = statMatrix(:,2);
% Enable the following to filter the data
[b1, a1] = butter(2, [59/500 61/500], 'stop');                      %% Remove 60hz Harmonic (noise)
goodDataTrace = filtfilt(b1, a1, goodDataTrace);                    %% Apply Filter
assignin('base', 'goodDataTrace', goodDataTrace);
badDataTrace = nan(size(goodDataTrace));
assignin('base', 'badDataTrace', badDataTrace);
end

%% Glenn notes
% Add a function that allows to switch channels
% May be better to use
% Have two vectors: original LFP and a NaN vector
% For all unwanted index ranges, swap values

% Look at Initialize plot script from Gabe
% AssignInBase
% For new files, make that the Y. data source