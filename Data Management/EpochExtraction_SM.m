function [unitEpoch, unitIDs, lfpEpoch, lfpIDs, trialTimeBins, eventTimeBins, trialInfo] = EpochExtraction_SM(eventRef, windowStart, windowEnd, varargin)
%% EpochExtraction_SM
% Runs through all the files in the current directory to extract user
% specified epochs from all the statMatrix files present and organize them
% into 3D structures for unit and LFP data as well as a vector containing
% event information.
%
%   Inputs:
%       'eventRef' : Reference point for the epoch extraction. Reference
%           events are limited to those used by the OrganizeTrialData_SM.m
%           function:
%               1) 'Odor' : Timepoint of odor presentation
%               2) 'PokeIn' : Timepoint triggering odor presentatino
%               3) 'PokeOut' : Timepoint of port withdrawal
%               4) 'FrontReward' : Timepoint of reward delivery 
%               5) 'RearReward' : Timepoint when reward was delivered at
%                   the back of the maze after a successfully completed
%                   sequence
%               6) 'ErrorSignal' : Timepoint when error signal happened.
%                   NOTE: there is no flag for this in the data collected
%                   in boston (initial CA1 data set) so don't use this flag
%                   for that data.
%       'windowStart' : Beginning of the epoch in seconds. This accepts 
%           negative values in order to extract a period prior to the event
%           reference point.
%       'windowEnd' : End point of the epoch in seconds.
%       'varargin' : Variable inputs to affect the way things run. Standard
%           structure of the input is a vector of cells with an identifier
%           string leading a parameter value.
%               Identifier:Values
%                   1) 'lfpBand': 'Raw'
%                               : 'Theta'
%                               : 'LowBeta'
%                               : 'Beta'
%                               : 'LowGamma'
%                               : 'HighGamma'
%                               : 'Ripple'
%                   2) 'lfpData': 'Phase' (the time series data is default,
%                                   'Phase' gives you the hilbert values);
%                               : 'Both' (This gives time series and phase)
%
%   Outputs:
%       'unitEpoch' : Extracted epoch with binary indicator for when a unit
%           is active. Organized in an MxNxO matrix where M = trial #, N =
%           timebin and O = unit.
%       'unitIDs' : Vector specifying the identity of the units stored in
%           the unitEpoch variable.
%       'lfpEpoch' : Extracted epoch of the LFP raw voltage data. Same
%           organization as the unitEpoch.
%       'lfpIDs' : Vector specifying the identity of source of the LFP 
%           signals.
%       'trialTimeBins' : Timebins for each trial.
%       'eventTimeBins' : Timebins relative to the extracted event.
%       'trialInfo' : Matrix containing information about the identity and
%           outcome of each trial. Rows represent individual trials;
%           columns are organized thusly:
%               Column 1) Trial performance. 1 = Correct, 0 = Incorrect
%               Column 2) Trial InSeq log. 1 = InSeq, 0 = OutSeq
%               Column 3) Trial Position.
%               Column 4) Trial Odor.
%
% Made by GE 08/14/2018

%% Locate files
files = dir(cd);
fileNames = {files.name};
matFiles = fileNames(cell2mat(cellfun(@(a)~isempty(a), strfind(fileNames, '.mat'), 'uniformoutput', 0)))';
if isempty(matFiles)
    matDir = uigetdir(cd, 'Select the folder with the statMatrix Files');
    cd(matDir);
end
statMatrixLog = false(size(matFiles));
for fl = 1:length(matFiles)
    variableInfo = who('-file', matFiles{fl});
    if sum(ismember(variableInfo, 'statMatrix'))==1
        statMatrixLog(fl) = true;
    end
end
smFiles = matFiles(statMatrixLog);

%% Load the behavior matrix file and create the eventAlignedMatrix
load(matFiles{cell2mat(cellfun(@(a)~isempty(a), strfind(matFiles, '_BehaviorMatrix'), 'uniformoutput', 0))});
eventAlignedMatrix = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [windowStart, windowEnd], eventRef);
trialInfo = [[eventAlignedMatrix.Performance];...
    [eventAlignedMatrix.TranspositionDistance]==0;...
    [eventAlignedMatrix.Position];...
    [eventAlignedMatrix.Odor]]';

%% Extract event timebins
eventTimes = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [0 0], eventRef);
trialTimeBins = ExtractTrialData_SM(eventAlignedMatrix, behavMatrix(:,1));
eventTimes = ExtractTrialData_SM(eventTimes, behavMatrix(:,1));
eventTimeBins = cellfun(@(a,b) a-b, trialTimeBins, eventTimes, 'uniformoutput',0);
frstNonMptTrl = find(cellfun(@(a)~isempty(a), eventTimeBins),1, 'first');
eventTimeBins = eventTimeBins{frstNonMptTrl};

fprintf('Trial and event timebins extracted\n');

%% Extract data from the statMatrix files
% Initialize the unit and LFP matrices
unitEpochs = cell(1,1,length(smFiles));
lfpEpochs = cell(1,1,length(smFiles));
unitIDs = [];
lfpIDs = [];

% Step through each file and extract the relevant data
for fl = 1:length(smFiles)
    load(smFiles{fl});
    fprintf('%s Loaded....', smFiles{fl});
    unitCols = cellfun(@(a)~isempty(a), regexp(statMatrixColIDs, 'T([0-9]*)-U([0-9]*)'));
    unitIDs = [unitIDs, statMatrixColIDs(unitCols)]; %#ok<AGROW>
    lfpCol = cellfun(@(a)~isempty(a), regexp(statMatrixColIDs, '_LFP_Raw\>'));
    lfpIDs = [lfpIDs, statMatrixColIDs(lfpCol)]; %#ok<AGROW>
    
    lfpEventData = ExtractTrialData_SM(eventAlignedMatrix, statMatrix(:,lfpCol)); %#ok<NODEF>
    
    mptLogLFP = cellfun(@(a)isempty(a), lfpEventData);
    lfpEventData(mptLogLFP) = {nan(size(eventTimeBins))};
    
    lfpEpochs{fl} = cell2mat(lfpEventData)';
    fprintf('LFP data collected......');
    numUnis = sum(unitCols);    
    fprintf('%i units found.....', numUnis);
    if numUnis >= 1
        uniEventData = ExtractTrialData_SM(eventAlignedMatrix,  statMatrix(:,unitCols));
        mptLogUni = cellfun(@(a)isempty(a), uniEventData);
        uniEventData(mptLogUni) = {nan(size(eventTimeBins,1), sum(unitCols))};
        tempUniData = nan(size(lfpEpochs{fl},1), size(lfpEpochs{fl},2), numUnis);
        for uni = 1:numUnis
            tempUniData(:,:,uni) = cell2mat(cellfun(@(a)a(:,uni), uniEventData, 'uniformoutput', 0))';
        end
        unitEpochs{fl} = tempUniData;
        fprintf('Unit data collected... File complete\n');
    else
        fprintf('File complete\n');
    end
end
%%
unitEpoch = cell2mat(unitEpochs);
lfpEpoch = cell2mat(lfpEpochs);

fprintf('Epoch extraction complete\n');
