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
trialWindow = [0, 0.5];

%% Evaluate Parameters
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
pokeInMatrix = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, trialWindow, 'PokeIn');
sampRate = 1/mode(diff(behavMatrix(:,1)));
corrWindow = sampRate./wSize;

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

%% Examine within band relations
for bnd = 1:length(freqBand)
    curBand = freqBand{bnd};
    curBandPower = cell2mat(reshape(powerVals(bnd,:), [1,1,size(powerVals,2)]));
    curBandPhase = cell2mat(reshape(phaseVals(bnd,:), [1,1,size(powerVals,2)]));
    regPowerCoupleTrl = nan(size(curBandPower,2),size(curBandPower,2),size(powerVals,2));
    regPhaseCohereTrl = nan(size(curBandPhase,2),size(curBandPhase,2),size(powerVals,2));
    for trl = 1:size(curBandPower,3)
        regPowerCoupleTrl(:,:,trl) = corr(curBandPower(:,:,trl));
        tempPhaseCohere = nan(size(curBandPhase,2));
        for reg1 = 1:size(curBandPhase,2)
            for reg2 = 1:size(curBandPhase,2)
                [tempPhaseCohere(reg1,reg2), ~] = circ_corrcc(curBandPhase(:,reg1,trl), curBandPhase(:,reg2,trl));
            end
        end
        regPhaseCohereTrl(:,:,trl) = tempPhaseCohere;
    end
    figure;
    subplot(2,2,1)
    imagesc(mean(regPowerCoupleTrl,3), [-1 1]);
    set(gca, 'ytick', 1:length(regions), 'yticklabel', regions,...
        'xtick', 1:length(regions), 'xticklabel', regions,...
        'DataAspectRatio',[1 1 1],'Layer','top');
    xtickangle(90);
    title([curBand ' Mean Regional Power Coupling ' num2str(trialWindow(1)) ' ' num2str(trialWindow(2))]);
    
    subplot(2,2,2)
    imagesc(median(regPowerCoupleTrl,3), [-1 1]);
    set(gca, 'ytick', 1:length(regions), 'yticklabel', regions,...
        'xtick', 1:length(regions), 'xticklabel', regions,...
        'DataAspectRatio',[1 1 1],'Layer','top');
    xtickangle(90);
    title([curBand ' Median Regional Power Coupling ' num2str(trialWindow(1)) ' ' num2str(trialWindow(2))]);
    
    subplot(2,2,3)
    imagesc(mean(regPhaseCohereTrl,3), [-1 1]);
    set(gca, 'ytick', 1:length(regions), 'yticklabel', regions,...
        'xtick', 1:length(regions), 'xticklabel', regions,...
        'DataAspectRatio',[1 1 1],'Layer','top');
    xtickangle(90);
    title([curBand ' Mean Regional Phase Coherence ' num2str(trialWindow(1)) ' ' num2str(trialWindow(2))]);
    
    subplot(2,2,4)
    imagesc(median(regPhaseCohereTrl,3), [-1 1]);
    set(gca, 'ytick', 1:length(regions), 'yticklabel', regions,...
        'xtick', 1:length(regions), 'xticklabel', regions,...
        'DataAspectRatio',[1 1 1],'Layer','top');
    xtickangle(90);
    title([curBand ' Median Regional Phase Coherence ' num2str(trialWindow(1)) ' ' num2str(trialWindow(2))]);
end
    



