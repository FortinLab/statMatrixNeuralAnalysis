function UnitSummary_SM
%% UnitSummary_SM
%
%   Omnibus code to summarize individual unit data

%%
close all
clear all

%% Determine Runtime Variables
pehBinSize = 0.05;
eventWindow = [-1 1];
spectFreqWindow = [1 120];
%%
origDir = cd;
[fileDir] = uigetdir(origDir);
if fileDir==0
    disp('Analysis Cancelled')
    return
else
    cd(fileDir)
end
dirContents = dir(fileDir);
fileNames = {dirContents.name};

%% Load Relevant Data
load(fileNames{cellfun(@(a)~isempty(a), strfind(fileNames, 'BehaviorMatrix'))});
load(fileNames{cellfun(@(a)~isempty(a), strfind(fileNames, 'EnsembleMatrix'))});

% Pull out the unitNames from the ensembleUnitSummaries variable.
% ***NOTE: if this variable does not exist the statMatrix data needs to be
% recompiled***
% We could pull this info from the ensembleMatrixColIDs, but doing it this
% way ensures that old data will be recompiled when run through this
% analysis.

% First pull out Unit names
unitNames = {ensembleUnitSummaries.UnitName};
% Then extract the prefix from the unit names which correspond to the
% tetrode the unit came from.
[tetStart, tetEnd] = regexp(unitNames, 'T([0-9]*)');
tetsWithUnits = unique(cellfun(@(a,b,c)a(b:c), unitNames, tetStart, tetEnd, 'uniformoutput', 0));

%% Create Behavior Matrices
pokeInAlignedBehavMatrix = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, eventWindow, 'PokeIn');
pokeOutAlignedBehavMatrix = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, eventWindow, 'PokeOut');
rewardAlignedBehavMatrix = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, eventWindow, 'FrontReward');
errorAlignedBehavMatrix = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, eventWindow, 'ErrorSignal');

behEventData = [pokeInAlignedBehavMatrix;...
    pokeOutAlignedBehavMatrix;...
    rewardAlignedBehavMatrix;...
    errorAlignedBehavMatrix];
behEventDataIDs = [{'Poke In'}, {'Poke Out'}, {'Reward'}, {'Error'}];

%% Create Trial Based Logical Vectors
seqLength = pokeInAlignedBehavMatrix(1).SeqLength;
allTrlLog = true(1,length(pokeInAlignedBehavMatrix));

% Performance vectors
corrTrlLog = [pokeInAlignedBehavMatrix.Performance]==1;
inCorrTrlLog = [pokeInAlignedBehavMatrix.Performance]==0;
perfLogs = [allTrlLog;corrTrlLog;inCorrTrlLog];
perfLogIDs = [{'All Trials'},{'Correct Trials'},{'Incorrect Trials'}];

% Temporal Context vectors
inSeqLog = [pokeInAlignedBehavMatrix.TranspositionDistance]==0;
outSeqLog = [pokeInAlignedBehavMatrix.TranspositionDistance]==1;
tcLogs = [allTrlLog; inSeqLog; outSeqLog];
tcLogIDs = [{'All Trials'}, {'InSeq Trials'}, {'OutSeq Trials'}];

% Odor vectors
for op = 1:seqLength
    odorVects.(sprintf('Odor%i', op)) = [pokeInAlignedBehavMatrix.Odor]==op;
    positionVects.(sprintf('Position%i', op)) = [pokeInAlignedBehavMatrix.Position]==op;
end

%% Unit Summary Overall
% Plot stuff here
for t = 1:length(tetsWithUnits)
    %% Find the relevant statMatrix data file in order to grab the relevant
    % LFP data
    curTet = tetsWithUnits{t};
    smFile = fileNames{cellfun(@(a)~isempty(a), strfind(fileNames, curTet))};
    load(smFile);
    
    curUnis = statMatrixColIDs(cellfun(@(a)~isempty(a), regexp(statMatrixColIDs, '-U([0-9]*)')));
    unitSummaryUnis = {unitSummary.UnitName};
    
    %% Plot Per-Tetrode Spectrograms
    % Performance
    PlotTrialEventSpect_SM(curTet, behEventData, behEventDataIDs,...
        perfLogs, perfLogIDs,...
        statMatrix, statMatrixColIDs, eventWindow, spectFreqWindow,0); %#ok<NODEF>
    
    % Temporal Context
    PlotTrialEventSpect_SM(curTet, behEventData, behEventDataIDs,...
        [allTrlLog;inSeqLog;outSeqLog], [{'All Trials'},{'InSeq Trials'},{'OutSeq Trials'}],...
        statMatrix, statMatrixColIDs, eventWindow, spectFreqWindow,0);
    
    %% Now do Per Unit Analysis
    for uni = 1:length(curUnis)
        curUnit = curUnis{uni};
        curUnitSummary = unitSummary(strcmp(curUnit, unitSummaryUnis));
        curUniSpikeLog = statMatrix(:,strcmp(curUnit,statMatrixColIDs))>0;
        curUniSpikeTimes = statMatrix(curUniSpikeLog ,1);
        
        %% Overall Summary Figures
        % Figure 1: Overall Unit Summary
        fig1 = figure('Name', [curUnit ' Summary'], 'NumberTitle', 'off');
        PlotUnitFeatures_SM(curUnitSummary, curUniSpikeTimes, fig1,0)
        
        %% Performance PEHs
        % Overall
        PlotTrialEventPEH_SM(curUnit, behEventData, behEventDataIDs,...
            perfLogs, perfLogIDs,...
            curUniSpikeLog, eventWindow, pehBinSize,0);
        
        % Performance by Odor
        PlotTrialEventByTrialTypePEH_SM(curUnit, behEventData, behEventDataIDs,...
            odorVects, fieldnames(odorVects),...
            perfLogs, perfLogIDs,...
            curUniSpikeLog, eventWindow, pehBinSize,0);
        
        % Performance by Position
        PlotTrialEventByTrialTypePEH_SM(curUnit, behEventData, behEventDataIDs,...
            positionVects, fieldnames(positionVects),...
            perfLogs, perfLogIDs,...
            curUniSpikeLog, eventWindow, pehBinSize,0);
                
        %% Temporal Context PEH
        % Overall
        PlotTrialEventPEH_SM(curUnit,  behEventData, behEventDataIDs,...
            tcLogs, tcLogIDs,...
            curUniSpikeLog, eventWindow, pehBinSize,0);
        
        % Temporal Context by Odor
        PlotTrialEventByTrialTypePEH_SM(curUnit, behEventData, behEventDataIDs,...
            odorVects, fieldnames(odorVects),...
            tcLogs, tcLogIDs,...
            curUniSpikeLog, eventWindow, pehBinSize,0);
        
        % Temporal Context by Position
        PlotTrialEventByTrialTypePEH_SM(curUnit, behEventData, behEventDataIDs,...
            positionVects, fieldnames(positionVects),...
            tcLogs, tcLogIDs,...
            curUniSpikeLog, eventWindow, pehBinSize,0);
        
        %% Spike-Phase Relations
        % To be fleshed out between GE and JDL
        
        % Overall Spike-Phase by Event
        
        % Spike-Phase by Event & Performance
        
        % Spike-Phase by Event & Temporal Context
        
        %%
        
        %%
        close all
    end
end
%%