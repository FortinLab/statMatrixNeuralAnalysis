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
% trialWindow = [-0.6 -0.1]; % Use to examine pre-trial period


colors = [[0 0 1];... %red
    [1 0 0];...   %blue
    [0 0.5 0];... %green
    [0.4 0 0.6];... %purple
    [1 0.5 0];... %orange
    [0 0.5 0.6]]; %yellow
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

%% Examine Overall Lag
allTrialSim = 1-corr(cell2mat(trialPopVects(perfLog)));
trialNums = [behavMatrixTrialStruct.TrialNum];
trialNums = trialNums(perfLog);
trlDiff = nan(length(trialNums));
for trlNumA = 1:length(trialNums)
    for trlNumB = 1:length(trialNums)
        trlDiff(trlNumA, trlNumB) = trialNums(trlNumB)-trialNums(trlNumA);
    end
end
meanSimByLagALL = nan(1,40);
stdSimByLagAll = nan(1,40);
for lag = 1:40
    lagLog = trlDiff==lag;
    meanSimByLagALL(lag) = mean(allTrialSim(lagLog));
    stdSimByLagALL(lag) = std(allTrialSim(lagLog))/sqrt(sum(lagLog(:))-1);
end
meanSimByLagALL = smooth(meanSimByLagALL./meanSimByLagALL(1))';
figure;
PlotLineAndFilledError(1:length(meanSimByLagALL), smooth(meanSimByLagALL)', stdSimByLagALL, [0 0.5 0.6]);
hold on;
%% Examine Sequence Lag
meanSimByLagODOR = nan(seqLength,40);
stdSimByLagODOR = nan(seqLength,40);

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
    
    curOdorMeanSimByLag = nan(1,40);
    curOdorStdSimByLag = nan(1,40);
    curOdorSEMSimByLag = nan(1,40);
    for lag = 1:40
        lagLog = tempSeqDiff==lag;
        curOdorMeanSimByLag(lag) = mean(curOdorSimMtx(lagLog));
        curOdorStdSimByLag(lag) = std(curOdorSimMtx(lagLog));
        curOdorSEMSimByLag(lag) = std(curOdorSimMtx(lagLog))/sqrt(sum(lagLog(:))-1);
    end
    meanSimByLagODOR(itm,:) = curOdorMeanSimByLag;
    stdSimByLagODOR(itm,:) = curOdorStdSimByLag;
    PlotLineAndFilledError(1:40, smooth(curOdorMeanSimByLag./curOdorMeanSimByLag(1))', smooth(curOdorSEMSimByLag)', colors(itm,:));
end

grandMeanSimByLagODOR = smooth(mean(meanSimByLagODOR./meanSimByLagODOR(:,1)));
PlotLineAndFilledError(1:40, grandMeanSimByLagODOR', smooth(std(meanSimByLagODOR./meanSimByLagODOR(:,1))/sqrt(size(meanSimByLagODOR,1)))', [0 0 0]);

% legend('Mean Sim By Trial #', 'OdorA', 'OdorB', 'OdorC', 'OdorD', 'OdorE', 'Odor Mean');
% legend('Mean Sim By Trial #', 'OdorA', 'OdorB', 'OdorC', 'OdorD', 'Odor Mean');
xlabel Lag
ylabel('Ensemble Distance (Normalized to Lag=1')
set(gca, 'ylim', [0.75 3]);
line(1:40, ones(1,40), 'linestyle', ':', 'linewidth', 2);

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
        

