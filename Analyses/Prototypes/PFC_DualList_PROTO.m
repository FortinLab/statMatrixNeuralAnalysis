
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

%% Define Standard Variables
slideWindowSize = 100;
spatialBinSize = 5;
numPerms = 50;

% Create Gaussian
instFRgauss = gausswin(slideWindowSize);
instFRgauss = instFRgauss/(length(instFRgauss)*mode(diff(behavMatrix(:,1))));

%% Load Relevant Data
load(fileNames{cellfun(@(a)~isempty(a), strfind(fileNames, 'BehaviorMatrix'))});
% load(fileNames{cellfun(@(a)~isempty(a), strfind(fileNames, 'EnsembleMatrix'))});

%%
pokeInTrialMatrix = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [-0.2 0.5], 'PokeIn');
pokeInTSs = behavMatrix(pokeInTrialMatrix(1).TrialLogVect,1)-behavMatrix(pokeInTrialMatrix(1).PokeInIndex,1);
pokeOutTrialMatrix = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [-0.5 0.2], 'PokeOut');
pokeOutTSs = behavMatrix(pokeOutTrialMatrix(1).TrialLogVect,1)-behavMatrix(pokeOutTrialMatrix(1).PokeOutIndex,1);

% Create Trial Logical Vectors
corrTrlLog = [pokeInTrialMatrix.Performance];
isLog = [pokeInTrialMatrix.TranspositionDistance]==0;
pos1log = [pokeInTrialMatrix.Position]==1;
% trialLog = (corrTrlLog & isLog & ~pos1log);
% seqType = 'All InSeq Correct Except Position 1';
trialLog = (corrTrlLog & isLog);
seqType = 'All InSeq Correct Trials';
% trialLog = (corrTrlLog);
% seqType = 'All Correct Trials';
posIDsISC = [pokeInTrialMatrix(trialLog).Position];
odorIDsISC = [pokeInTrialMatrix(trialLog).Odor];
%%