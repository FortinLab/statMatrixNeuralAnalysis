function UnitSummary_SM
%% UnitSummary_SM
%
%   Omnibus code to summarize individual unit data

%%
close all
clear all

%% Determine Runtime Variables
pehBinSize = 0.1;
eventWindow = [-0.5 0.5];
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

%% Create Trial Based Logical Vectors
seqLength = pokeInAlignedBehavMatrix(1).SeqLength;

% Performance vectors
corrTrlLog = [pokeInAlignedBehavMatrix.Performance]==1;
inCorrTrlLog = [pokeInAlignedBehavMatrix.Performance]==0;

% Temporal Context vectors
inSeqLog = [pokeInAlignedBehavMatrix.TranspositionDistance]==0;
outSeqLog = [pokeInAlignedBehavMatrix.TranspositionDistance]==1;

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
    
    for uni = 1:length(curUnis)
        curUnit = curUnis{uni};
        curUnitSummary = unitSummary(strcmp(curUnit, unitSummaryUnis));
        curUniSpikeLog = statMatrix(:,strcmp(curUnit,statMatrixColIDs))>0; %#ok<NODEF>
        curUniSpikeTimes = statMatrix(curUniSpikeLog ,1);
        
        %% Overall Summary Figures
        % Figure 1: Overall Unit Summary
        fig1 = figure('Name', [curUnit ' Summary'], 'NumberTitle', 'off');
        PlotUnitFeatures_SM(curUnitSummary, curUniSpikeTimes, fig1)
        
        % Figure 2: Spectrogram & epoc power (GG & GE to include)
        % fig2 = figure; 
        
        %% Figure 3: PEH Overall Perf
        pehOvrlPerf = figure('Name', [curUnit ' Trial Event by Performance'], 'NumberTitle', 'off');
        PlotTrialEventPEH_SM(curUnit, pokeInAlignedBehavMatrix, pokeOutAlignedBehavMatrix, rewardAlignedBehavMatrix, errorAlignedBehavMatrix,...
            corrTrlLog, 'Correct Trials', inCorrTrlLog, 'Incorrect Trials',...
            curUniSpikeLog, eventWindow, pehBinSize, pehOvrlPerf);        
        
        %% Figure 4: PEH Temporal Context
        pehOvrlTC = figure('Name', [curUnit ' Trial Event by Temporal Context'], 'NumberTitle', 'off');
        PlotTrialEventPEH_SM(curUnit, pokeInAlignedBehavMatrix, pokeOutAlignedBehavMatrix, rewardAlignedBehavMatrix, errorAlignedBehavMatrix,...
            inSeqLog, 'InSeq Trials', outSeqLog, 'OutSeq Trials',...
            curUniSpikeLog, eventWindow, pehBinSize, pehOvrlTC);
        
        %% Figures: PEH Performance by Odor
        PlotTrialEventByTrialTypePEH_SM(curUnit, [pokeInAlignedBehavMatrix; pokeOutAlignedBehavMatrix; rewardAlignedBehavMatrix; errorAlignedBehavMatrix],...
            [{'Poke In'}, {'Poke Out'}, {'Reward'}, {'Error'}],...
            odorVects, fieldnames(odorVects),...
            [true(1,size(inSeqLog,2));corrTrlLog;inCorrTrlLog], [{'All Trials'},{'Correct Trials'},{'Incorrect Trials'}],...
            curUniSpikeLog, eventWindow, pehBinSize);
        
        %% Figures: PEH Temporal Context by Odor
        PlotTrialEventByTrialTypePEH_SM(curUnit, [pokeInAlignedBehavMatrix; pokeOutAlignedBehavMatrix; rewardAlignedBehavMatrix; errorAlignedBehavMatrix],...
            [{'Poke In'}, {'Poke Out'}, {'Reward'}, {'Error'}],...
            odorVects, fieldnames(odorVects),...
            [true(1,size(inSeqLog,2));inSeqLog;outSeqLog], [{'All Trials'},{'InSeq Trials'},{'OutSeq Trials'}],...
            curUniSpikeLog, eventWindow, pehBinSize);
        
        %% Figures: PEH Performance by Position
         PlotTrialEventByTrialTypePEH_SM(curUnit, [pokeInAlignedBehavMatrix; pokeOutAlignedBehavMatrix; rewardAlignedBehavMatrix; errorAlignedBehavMatrix],...
            [{'Poke In'}, {'Poke Out'}, {'Reward'}, {'Error'}],...
            positionVects, fieldnames(positionVects),...
            [true(1,size(inSeqLog,2));corrTrlLog;inCorrTrlLog], [{'All Trials'},{'Correct Trials'},{'Incorrect Trials'}],...
            curUniSpikeLog, eventWindow, pehBinSize);
                 
        %% Figures: PEH Temporal Context by Position
        PlotTrialEventByTrialTypePEH_SM(curUnit, [pokeInAlignedBehavMatrix; pokeOutAlignedBehavMatrix; rewardAlignedBehavMatrix; errorAlignedBehavMatrix],...
            [{'Poke In'}, {'Poke Out'}, {'Reward'}, {'Error'}],...
            positionVects, fieldnames(positionVects),...
            [true(1,size(inSeqLog,2));inSeqLog;outSeqLog], [{'All Trials'},{'InSeq Trials'},{'OutSeq Trials'}],...
            curUniSpikeLog, eventWindow, pehBinSize);
        %% Figures organized by position
        % PEH By Position Poke In Perf
%         piPos = figure;
        
        % PEH By Position Poke Out Perf
%         poPos = figure;
        
        % PEH By Position Reward Perf
%         rwdPos = figure;
        
        % PEH By Position Error Perf
%         errPos = figure
        
        
        
        %
        
        % Figure 5: PEH By Position Perf
%         fig5 = figure;
        
        % Figure 6: PEH Overall TC
%         fig6 = figure;
        
        % Figure 7: PEH By Odor TC
%         fig7 = figure;
        
        % Figure 8: PEH By Position TC
%         fig8 = figure;
        
        % Figure 9: Spike-Phase Overall Perf (JDL to include)
        % fig9 = figure;
        
        % Figure 10: Spike-Phase By Odor Perf (JDL)
        % fig10 = figure;
        
        % Figure 11: Spike-Phase By Position Perf (JDL)
        % fig11 = figure;
        
        % Figure 12: Spike-Phase Overall TC (JDL)
        % fig12 = figure;
        
        % Figure 13: Spike-Phase By Odor TC (JDL)
        % fig13 = figure;
        
        % Figure 14: Spike-Phase By Position TC (JDL)
        % fig14 = figure;
        
        % Figure 15: Something magical?.... there will probably be a 15
        % fig15 = figure;
        %%
        close all
    end
end
%%