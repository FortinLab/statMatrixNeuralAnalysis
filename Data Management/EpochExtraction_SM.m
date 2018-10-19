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
%               2) 'PokeIn' : Timepoint triggering odor presentation
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
%           string leading a parameter value. Defaults reflect what data is
%           passed were the identifiers not passed as varargin
%               *Identifier:Values*
%                   1) 'lfpBand': Default is Raw time series data
%                               : 'Theta' - 4:12Hz filtered
%                               : 'LowBeta' - 13:19Hz filtered
%                               : 'Beta' - 20:40Hz filtered
%                               : 'LowGamma' - 41:59Hz filtered
%                               : 'HighGamma' - 59:61Hz filtered
%                               : 'Ripple' - 150:250Hz filtered
%                               : 'All' - Pulls all the LFP band columns
%                   2) 'lfpData': Default = voltage values
%                               : 'Phase' - Hilbert values
%                               : 'Both' - voltage and Hilbert values
%                   3) 'org'    : Default = 'TrTiU'
%                               : 'TrTiU' - Organizes the unit/lfp epoc 
%                                   data in a Trial X Time X Unit matrix
%                               : 'TiUTr' - Organizes the unit/lfp epoc
%                                   data in a Time X Unit X Trial matrix
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

%% Evaluate varargin inputs
if sum(strcmp(varargin, 'lfpBand'))==1
    lfpBandSpot = find(strcmp(varargin, 'lfpBand'));
    bandIDs = varargin{lfpBandSpot+1};
    if iscell(bandIDs)
        bandID = '([';
        b = 1;
        while b<=length(bandIDs)
            bandID = [bandID '|' bandIDs{b}]; %#ok<AGROW>
            b = b+1;
        end
        bandID = [bandID '])'];
    elseif strcmp(bandIDs, 'All')
        bandID = '([A-Z | a-z]*)';
    else
        bandID = bandIDs;
    end
else
    bandID = 'Raw';    
end
if sum(strcmp(varargin, 'lfpData'))==1
    lfpDataSpot = find(strcmp(varargin, 'lfpData'));
    dataID = varargin{lfpDataSpot+1};
    if strcmp(dataID, 'Both')
        lfpColRegXprsnString = ['_LFP_' bandID];
    elseif strcmp(dataID, 'Phase')
        lfpColRegXprsnString = ['_LFP_' bandID '_HilbVals'];
    end
else
    lfpColRegXprsnString = ['_LFP_' bandID '\>'];
end
if sum(strcmp(varargin, 'org'))==1
    orgDataSpot = find(strcmp(varargin, 'org'));
    dataOrgID = varargin{orgDataSpot+1};
    if strcmp(dataOrgID, 'TrTiU')
        org = 1;
    elseif strcmp(dataOrgID, 'TiUTr')
        org = 2;
    end
else
    org = 1;
end
%% Extract data from the statMatrix files
% Initialize the unit and LFP matrices
if org == 1
    unitEpochs = cell(1,1,length(smFiles));
    lfpEpochs = cell(1,1,length(smFiles));
    unitIDs = [];
    lfpIDs = [];
    
    % Step through each file and extract the relevant data
    for fl = 1:length(smFiles)
        load(smFiles{fl});
        fprintf('%s Loaded....', smFiles{fl});
        unitCols = cellfun(@(a)~isempty(a), regexp(statMatrixColIDs, '-U([0-9]*)'));
        unitIDs = [unitIDs, statMatrixColIDs(unitCols)]; %#ok<AGROW>
        lfpCol = cellfun(@(a)~isempty(a), regexp(statMatrixColIDs, lfpColRegXprsnString));
        lfpIDs = [lfpIDs, statMatrixColIDs(lfpCol)]; %#ok<AGROW>
        
        lfpEventData = ExtractTrialData_SM(eventAlignedMatrix, statMatrix(:,lfpCol)); %#ok<NODEF>
        mptLogLFP = cellfun(@(a)isempty(a), lfpEventData);
        lfpEventData(mptLogLFP) = {nan(size(eventTimeBins))};
        lfpEventData = cellfun(@(a)reshape(a, [1, size(a,1), size(a,2)]), lfpEventData, 'uniformoutput',0);
        
        lfpEpochs{1,1,fl} = cell2mat(lfpEventData');
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
            unitEpochs{1,1,fl} = tempUniData;
            fprintf('Unit data collected... File complete\n');
        else
            fprintf('File complete\n');
        end
    end
elseif org == 2
    unitEpochs = cell(1,length(smFiles),length(eventAlignedMatrix));
    lfpEpochs = cell(1,length(smFiles),length(eventAlignedMatrix));
    unitIDs = [];
    lfpIDs = [];
    % Step through each file and extract the relevant data
    for fl = 1:length(smFiles)
        load(smFiles{fl});
        fprintf('%s Loaded....', smFiles{fl});
        unitCols = cellfun(@(a)~isempty(a), regexp(statMatrixColIDs, '-U([0-9]*)'));
        unitIDs = [unitIDs, statMatrixColIDs(unitCols)]; %#ok<AGROW>
        lfpCol = cellfun(@(a)~isempty(a), regexp(statMatrixColIDs, lfpColRegXprsnString));
        lfpIDs = [lfpIDs, statMatrixColIDs(lfpCol)]; %#ok<AGROW>
        
        lfpEventData = ExtractTrialData_SM(eventAlignedMatrix, statMatrix(:,lfpCol)); %#ok<NODEF>
        mptLogLFP = cellfun(@(a)isempty(a), lfpEventData);
        lfpEventData(mptLogLFP) = {nan(size(eventTimeBins))};
        lfpEpochs(:,fl,:) = lfpEventData;
        fprintf('LFP data collected......');
        
        numUnis = sum(unitCols);
        fprintf('%i units found.....', numUnis);
        if numUnis >= 1
            uniEventData = ExtractTrialData_SM(eventAlignedMatrix,  statMatrix(:,unitCols));
            mptLogUni = cellfun(@(a)isempty(a), uniEventData);
            uniEventData(mptLogUni) = {nan(size(eventTimeBins,1), sum(unitCols))};            
            unitEpochs(:,fl,:) = uniEventData;
            fprintf('Unit data collected... File complete\n');
        else
            fprintf('File complete\n');
        end
    end
    unitEpochs(:,cellfun(@(a)isempty(a),unitEpochs(1,:,1)),:) = [];
end
%%
unitEpoch = cell2mat(unitEpochs);
lfpEpoch = cell2mat(lfpEpochs);

fprintf('Epoch extraction complete\n');
