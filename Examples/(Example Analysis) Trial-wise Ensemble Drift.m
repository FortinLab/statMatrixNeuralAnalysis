%% Trial-wise Ensemble Drift Analysis %%
% This analysis is included to work as an example to demonstrate how the
% statMatrix data created for a recording session could be used in an
% analysis. 
% 
% In this analysis I'm examining the population vectors derived from
% different trial periods and comparing them as a function of their
% distance from one another during the session to see how similar those
% population vectors are.
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

%%
load(fileNames{cellfun(@(a)~isempty(a), strfind(fileNames, 'BehaviorMatrix'))});
load(fileNames{cellfun(@(a)~isempty(a), strfind(fileNames, 'EnsembleMatrix'))});

%% Runtime Variables
trialWindow = [0 1]; % Use to examine trial period
% trialWindow = [-0.5 0]; % Use to examine pre-trial period


colors = [[110/255 190/255 245/255];... %A
    [180/255 170/255 160/255];...   %B
    [69/255 194/255 135/255];... %C
    [164/255 120/255 189/255];... %D
    [254/255 149/255 95/255]]; %E
%% Organize Data by Trials & Extract Relevant Data
% Create Behavior structure
behavMatrixTrialStruct = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, trialWindow, 'PokeIn');
ensembleTrialOrg = ExtractTrialData_SM(behavMatrixTrialStruct, ensembleMatrix(:,2:end));
trialPopVects = cellfun(@(a){sum(a)'}, ensembleTrialOrg);

seqLength = behavMatrixTrialStruct(1).SeqLength;
seqLength = 4;

% Extract Trial Odor IDs
trialOdor = [behavMatrixTrialStruct.Odor];
% Extract Trial Positions
trialPosition = [behavMatrixTrialStruct.Position];
% Extract Performance Log
perfLog = [behavMatrixTrialStruct.Performance]==1;
% Extract Sequence Number
sequenceNum = [behavMatrixTrialStruct.SequenceNum];
% Identify InSeq Trials
inSeqLog = trialOdor==trialPosition;

%% Examine Overall Lag
trlLog = perfLog & inSeqLog;
% trlLog = perfLog;
allTrialSim = 1-corr(cell2mat(trialPopVects(trlLog)));
trialNums = [behavMatrixTrialStruct.TrialNum];
trialNums = trialNums(trlLog);
trlDiff = nan(length(trialNums));
for trlNumA = 1:length(trialNums)
    for trlNumB = 1:length(trialNums)
        trlDiff(trlNumA, trlNumB) = trialNums(trlNumB)-trialNums(trlNumA);
    end
end
meanSimByLagALL = nan(1,30);
stdSimByLagALL = nan(1,30);
lagPos = 0;
for lag = 1:seqLength:(30*seqLength)
    lagPos = lagPos+1;
    lagLog = trlDiff==lag;
    meanSimByLagALL(lagPos) = mean(allTrialSim(lagLog));
    stdSimByLagALL(lagPos) = std(allTrialSim(lagLog))/sqrt(sum(lagLog(:))-1);
end
meanSimByLagALL = smooth(meanSimByLagALL./meanSimByLagALL(1))';
figure;
PlotLineAndFilledError(1:length(meanSimByLagALL), meanSimByLagALL, stdSimByLagALL, [0 0 0]);
hold on;
%% Examine Sequence Lag
meanSimByLagODOR = nan(seqLength,30);
stdSimByLagODOR = nan(seqLength,30);

for itm = 1:seqLength
    odorLog = trialOdor==itm;
    posLog = trialPosition==itm;
    trialLog = odorLog & posLog & perfLog; 
    
    curPopVects = cell2mat(trialPopVects(trialLog));
    curOdorSimMtx = 1-corr(curPopVects);
    curSeqNums = [behavMatrixTrialStruct(trialLog).SequenceNum];
    tempSeqDiff = nan(length(curSeqNums));
    for seqNumA = 1:length(curSeqNums)
        for seqNumB = 1:length(curSeqNums)
            tempSeqDiff(seqNumA, seqNumB) = curSeqNums(seqNumB)-curSeqNums(seqNumA);
        end
    end
    
    curOdorMeanSimByLag = nan(1,30);
    curOdorStdSimByLag = nan(1,30);
    curOdorSEMSimByLag = nan(1,30);
    for lag = 1:30
        lagLog = tempSeqDiff==lag;
        curOdorMeanSimByLag(lag) = mean(curOdorSimMtx(lagLog));
        curOdorStdSimByLag(lag) = std(curOdorSimMtx(lagLog));
        curOdorSEMSimByLag(lag) = std(curOdorSimMtx(lagLog))/sqrt(sum(lagLog(:))-1);
    end
    meanSimByLagODOR(itm,:) = curOdorMeanSimByLag;
    stdSimByLagODOR(itm,:) = curOdorStdSimByLag;
    PlotLineAndFilledError(1:30, smooth(curOdorMeanSimByLag./curOdorMeanSimByLag(1))', smooth(curOdorSEMSimByLag)', colors(itm,:));
end

grandMeanSimByLagODOR = smooth(mean(meanSimByLagODOR./meanSimByLagODOR(:,1)));
plot(1:30, grandMeanSimByLagODOR', 'color', 'black', 'linewidth', 2, 'linestyle', '--');

xlabel Lag
ylabel('Ensemble Distance (Normalized to Lag=1')
set(gca, 'ylim', [0.75 3]);
line(1:30, ones(1,30), 'linestyle', ':', 'linewidth', 2);

%% Examine Lag Effects Around OutSeq Trials
inSeqLog = trialOdor == trialPosition;

% Identify OutSeq Correct Trials
outSeqCorrTrls = find(~inSeqLog & perfLog);

% Create Data Storage Structure
driftRatioA = nan(1,length(outSeqCorrTrls));
driftAzScore = nan(1,length(outSeqCorrTrls));

for osTrl = 1:length(outSeqCorrTrls)
    curTrl = outSeqCorrTrls(osTrl);
    curTrlPos = trialPosition(curTrl);
    curSeqNum = sequenceNum(curTrl);
    prevSeqTrlPos = find(sequenceNum==(curSeqNum-1) & trialPosition==curTrlPos);
    if ~isempty(prevSeqTrlPos) &... % IF the rat reached the same trial on the previous sequence
            (curTrlPos < seqLength) &... % AND the outSeq trial was not the last in the sequence
            perfLog(curTrl-1)==1 &... % AND the inSeq trial in the position before it was correct (it should be in almost all cases)
            perfLog(prevSeqTrlPos-1)==1 &... % AND the inSeq trial in the position before it from the previous sequence was correct
            perfLog(curTrl+1)==1 &... % AND the inSeq trial in the position following it was correct
            perfLog(prevSeqTrlPos+1)==1 % AND the inSeq trial in the position following from the previous sequence was correct
        
        prevTrlNorm = meanSimByLagODOR(curTrlPos-1,1);
        prevSeqPrevTrlEnsemble = trialPopVects{prevSeqTrlPos-1};
        curSeqPrevTrlEnsemble = trialPopVects{curTrl-1};
        prevTrlCorr = 1-corr(prevSeqPrevTrlEnsemble, curSeqPrevTrlEnsemble);
        
        subseqTrlNorm = meanSimByLagODOR(curTrlPos+1,1);
        prevSeqSubseqTrlEnsemble = trialPopVects{prevSeqTrlPos+1};
        curSeqSubseqTrlEnsemble = trialPopVects{curTrl+1};
        subseqTrlCorr = 1-corr(prevSeqSubseqTrlEnsemble, curSeqSubseqTrlEnsemble);
        
        driftRatioA(osTrl) = subseqTrlCorr/prevTrlCorr;
        driftAzScore(osTrl) = (subseqTrlCorr-meanSimByLagODOR(curTrlPos+1,1))/stdSimByLagODOR(curTrlPos+1,1);
%         driftRatioA(osTrl) = (subseqTrlCorr/prevTrlCorr)/(subseqTrlNorm/prevTrlNorm);
    end
end
figure;
boxplot(driftRatioA(~isnan(driftRatioA)))
hold on; 
scatter(ones(1,sum(~isnan(driftRatioA))), driftRatioA(~isnan(driftRatioA)), 'black*');
title('OutSeq Influence on Drift')
set(gca, 'xticklabel', []);
ylabel([{'Drift Ratio'}; {'(Trial after OutSeq / Trial before OutSeq)'}]);
        

