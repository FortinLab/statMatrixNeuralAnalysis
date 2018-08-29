function MultiSiteNetworkMeasures_PROTO
%% MultiSiteNetworkMeasures_PROTO
%   Prototype for network measures based on inter-regional correlations of
%   LFP signals.
%   Analyses focused on examining phase and power coupling between regions
%   to identify dynamics of information flow
%
%% Define Analysis Parameters
% Theta, LowBeta, Beta, LowGamma, HighGamma, Ripple
minFreqs = num2cell([4,13,20,41,61,150]);
maxFreqs = num2cell([12,19,40,59,80,250]);
freqBand = [{'Theta'}, {'LowBeta'}, {'Beta'}, {'LowGamma'}, {'HighGamma'}, {'Ripple'}];
wSize = cellfun(@(a,b) 1.25*mean([a,b]), minFreqs, maxFreqs);

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

%% Load the behavior matrix and create organizational structures
load(fileNames{cellfun(@(a)~isempty(a), strfind(fileNames, 'BehaviorMatrix'))});
pokeInMatrix = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [-0.5, 1.5], 'PokeIn');

%% Extract Power and Phase Values
smFiles = matFiles(statMatrixLog);
file1info = strsplit(smFiles{1},'_');
fileInfo = cell(length(smFiles), length(file1info));
powerVals = cell(length(freqBand), length(pokeInMatrix));
phaseVals = cell(length(freqBand), length(pokeInMatrix));
for fl = 1:length(smFiles)
    fileInfo(fl,:) = strsplit(smFiles{fl}, '_');
    fprintf('Running %s...', smFiles{fl});
    load(smFiles{fl});
    for bnd = 1:length(freqBand)
        curBandData = statMatrix(:, cellfun(@(a)~isempty(a), regexp(statMatrixColIDs, ['_LFP_' freqBand{bnd} '\>']))); %#ok<NODEF>
        [tempPowerVals,~] = mvrms_mfile_simulink(statMatrix(:,1),curBandData,wSize(bnd),0);
        tempPowerValsByTrial = ExtractTrialData_SM(pokeInMatrix, zscore(tempPowerVals));
        powerVals(bnd,:) = cellfun(@(a,b)horzcat(a,b), powerVals(bnd,:), tempPowerValsByTrial, 'uniformoutput', 0);

        tempPhaseVals = statMatrix(:,cellfun(@(a)~isempty(a), regexp(statMatrixColIDs, ['_LFP_' freqBand{bnd} '_HilbVals\>'])));
        tempPhaseValsByTrial = ExtractTrialData_SM(pokeInMatrix, tempPhaseVals);
        phaseVals(bnd,:) = cellfun(@(a,b)horzcat(a,b), phaseVals(bnd,:), tempPhaseValsByTrial, 'uniformoutput', 0);
        fprintf('%s done...', freqBand{bnd});
    end
    fprintf('\n');
end
regions = fileInfo(:,end-2);
chans = fileInfo(:,end-1);



