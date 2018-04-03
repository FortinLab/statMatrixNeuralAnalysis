function UnitSummary_SM
%% UnitSummary_SM
%
%   Omnibus code to summarize individual unit data

%%
close all
clear all
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

%% Create Behavior Matrices
pokeInAlignedBehavMatrix = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [-0.5 0.5], 'PokeIn');
ensemblePokeInOrg = ExtractTrialData_SM(pokeInAlignedBehavMatrix, ensembleMatrix(:,2:end));

pokeOutAlignedBehavMatrix = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [-0.5 0.5], 'PokeOut');
ensemblePokeOutOrg = ExtractTrialData_SM(pokeOutAlignedBehavMatrix, ensembleMatrix(:,2:end));

rewardAlignedBehavMatrix = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [-0.5 0.5], 'FrontReward');
ensembleRewardOrg = ExtractTrialData_SM(rewardAlignedBehavMatrix, ensembleMatrix(:,2:end));

errorAlignedBehavMatrix = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [-0.5 0.5], 'ErrorSignal');
ensembleErrorOrg = ExtractTrialData_SM(errorAlignedBehavMatrix, ensembleMatrix(:,2:end));


%% Unit Summary Overall
% Plot stuff here

%%