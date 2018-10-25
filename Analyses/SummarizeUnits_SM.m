function SummarizeUnits_SM
%%
%%
files = dir(cd);
fileNames = {files.name};
matFiles = fileNames(cell2mat(cellfun(@(a)~isempty(a), strfind(fileNames, '.mat'), 'uniformoutput', 0)))';
if isempty(matFiles)
    matDir = uigetdir(cd, 'Select the folder with the statMatrix Files');
    if matDir==0
        disp('Analysis Cancelled')
        return
    else
        cd(matDir)
    end
end
statMatrixLog = false(size(matFiles));
for fl = 1:length(matFiles)
    variableInfo = who('-file', matFiles{fl});
    if sum(ismember(variableInfo, 'statMatrix'))==1
        statMatrixLog(fl) = true;
    end
end
smFiles = matFiles(statMatrixLog);
dirContents = dir(cd);
fileNames = {dirContents.name};

%% Load Ensemble and Behavior Data Matrices
load(fileNames{cellfun(@(a)~isempty(a), strfind(fileNames, 'BehaviorMatrix'))});
load(fileNames{cellfun(@(a)~isempty(a), strfind(fileNames, 'EnsembleMatrix'))});

%% Create data storage structures and standard variables
unitIDs = ensembleMatrixColIDs(2:end); %#ok<COLND>

unitInfo = ensembleUnitSummaries;
for u = 1:length(unitInfo)
    unitInfo(u).Directory = cd;
end

% Window size for sliding window analyses
slideWindowSize = 200;

% Gaussian used for calculating instantaneous firing rate
instFRgauss = gausswin(slideWindowSize);
instFRgauss = instFRgauss/(length(instFRgauss)*mode(diff(behavMatrix(:,1))));

% Phase bins used for LFP phase analysis
phaseBins = linspace(-pi,pi,4);
phaseBinLabels = cell(1,length(phaseBins)-1);
for p = 1:length(phaseBins)-1
    if phaseBins(p)<0
        lowBin = sprintf('neg%1.0g', abs(phaseBins(p)));
    else
        lowBin = sprintf('%1.0g', phaseBins(p));
    end
    if phaseBins(p+1)<0
        highBin = sprintf('neg%1.0g', abs(phaseBins(p+1)));
    else
        highBin = sprintf('%1.0g', phaseBins(p+1));
    end
    phaseBinLabels{p} = sprintf('From_%s_to_%s', lowBin, highBin);
end

% Number of permutations for chance distribution estimation
numPerms = 100;

%% Analysis #1: Examine Average Evoked activity during trial periods
fprintf('Starting Analysis #1....');
preTrialBehavMatrix = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [-0.5 0], 'PokeIn');
preTrialEnsemble = ExtractTrialData_SM(preTrialBehavMatrix, ensembleMatrix(:,2:end)); %#ok<*NODEF>

earlyTrialBehavMatrix = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [0 0.5], 'PokeIn');
earlyTrialEnsemble = ExtractTrialData_SM(earlyTrialBehavMatrix, ensembleMatrix(:,2:end)); 

lateTrialBehavMatrix = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [-0.5 0], 'PokeOut');
lateTrialEnsemble = ExtractTrialData_SM(lateTrialBehavMatrix, ensembleMatrix(:,2:end)); 

postTrialBehavMatrix = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [0 0.5], 'PokeOut');
postTrialEnsemble = ExtractTrialData_SM(postTrialBehavMatrix, ensembleMatrix(:,2:end)); 

rewardBehavMatrix = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [-0.25 0.25], 'FrontReward');
rewardEnsemble = ExtractTrialData_SM(rewardBehavMatrix, ensembleMatrix(:,2:end)); 

errorBehavMatrix = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [0 0.5], 'ErrorSignal');
errorEnsemble = ExtractTrialData_SM(errorBehavMatrix, ensembleMatrix(:,2:end)); 

% Extract Trial Information
correctTrialLog = logical([preTrialBehavMatrix.Performance]);
transDist = [preTrialBehavMatrix.TranspositionDistance];
inSeqTrialLog = transDist==0;
trialOdor = [preTrialBehavMatrix.Odor];
trialPosition = [preTrialBehavMatrix.Position];
% Create more specialized logical vectors 
transInLog = false(size(correctTrialLog));
prevOdrLog = false(size(correctTrialLog));
prevOdrVect = nan(size(correctTrialLog));
curPosVect = nan(size(correctTrialLog));
nextPosLog = false(size(correctTrialLog));
nextPosVect = nan(size(correctTrialLog));
curOdrVect = nan(size(correctTrialLog));
for trl = 2:length(transInLog)
    if ~inSeqTrialLog(trl-1) && trialPosition(trl)-trialPosition(trl-1)==1
        transInLog(trl) = true;
    end
    if correctTrialLog(trl) 
        if trialPosition(trl)-trialPosition(trl-1)==1
            prevOdrLog(trl) = true;
            prevOdrVect(trl) = trialOdor(trl-1);
            curPosVect(trl) = trialPosition(trl);
        end
        if trl+1<length(transInLog) && correctTrialLog(trl+1) && trialPosition(trl+1)-trialPosition(trl)==1
            nextPosLog(trl) = true;
            nextPosVect(trl) = trialPosition(trl+1);
            curOdrVect(trl) = trialOdor(trl);
        end
    end
end


for u = 1:length(unitIDs)
    unitInfo(u).TrialInfo.CorrectTrialLog = correctTrialLog;
    unitInfo(u).TrialInfo.InSeqTrialLog = inSeqTrialLog;
    unitInfo(u).TrialInfo.Odor = trialOdor;
    unitInfo(u).TrialInfo.Position = trialPosition;
    unitInfo(u).TrialInfo.TransInLog = transInLog;
    unitInfo(u).TrialInfo.TransDist = transDist;
    unitInfo(u).TrialInfo.PokeDuration = [preTrialBehavMatrix.PokeDuration];
    
    % Extract Unit Firing Activity
    curUniPreTrialActivity = cell2mat(cellfun(@(a)sum(a(:,u)), preTrialEnsemble, 'uniformoutput', 0));    
    curUniErlyTrialActivity = cell2mat(cellfun(@(a)sum(a(:,u)), earlyTrialEnsemble, 'uniformoutput', 0));        
    curUniLtTrialActivity = cell2mat(cellfun(@(a)sum(a(:,u)), lateTrialEnsemble, 'uniformoutput', 0));    
    curUniPstTrlActivity = cell2mat(cellfun(@(a)sum(a(:,u)), postTrialEnsemble, 'uniformoutput', 0));    
    curUniRwdTrlActivity = cell2mat(cellfun(@(a)sum(a(:,u)), rewardEnsemble, 'uniformoutput', 0));    
    curUniErrTrlActivity = cell2mat(cellfun(@(a)sum(a(:,u)), errorEnsemble, 'uniformoutput', 0));
    
    unitInfo(u).TrialEpochFRs.PreTrialFR = [mean(curUniPreTrialActivity) std(curUniPreTrialActivity)];
    unitInfo(u).TrialEpochFRs.EarlyTrialFR = [mean(curUniErlyTrialActivity) std(curUniErlyTrialActivity)];
    unitInfo(u).TrialEpochFRs.LateTrialFR = [mean(curUniLtTrialActivity) std(curUniLtTrialActivity)];
    unitInfo(u).TrialEpochFRs.PostTrialFR = [mean(curUniPstTrlActivity) std(curUniPstTrlActivity)];
    unitInfo(u).TrialEpochFRs.RewardFR = [mean(curUniRwdTrlActivity(correctTrialLog)) std(curUniRwdTrlActivity(correctTrialLog))];
    unitInfo(u).TrialEpochFRs.ErrorFR = [mean(curUniErrTrlActivity(~correctTrialLog)) std(curUniErrTrlActivity(~correctTrialLog))];
    
    % Evaluate Unit Firing Across Time Periods
    frPeriodDta = [curUniPreTrialActivity(correctTrialLog), curUniErlyTrialActivity(correctTrialLog), curUniLtTrialActivity(correctTrialLog), curUniPstTrlActivity(correctTrialLog)]';
    frPeriodGrp = [ones(1,sum(correctTrialLog)), ones(1,sum(correctTrialLog))*2, ones(1,sum(correctTrialLog))*3, ones(1,sum(correctTrialLog))*4]';
    [~,frPeriodTable,~] = anova1(frPeriodDta, frPeriodGrp, 'off');
    unitInfo(u).TrialEpochStats.TrialEpochsF = [frPeriodTable{2,5}, frPeriodTable{2,6}];
    
    % Evaluate Unit Feeback Properties
    feedbackDta = [curUniRwdTrlActivity(correctTrialLog), curUniErrTrlActivity(~correctTrialLog)];
    feedbackGrps = [ones(1,sum(correctTrialLog)), ones(1,sum(~correctTrialLog))*2];
    [~,feedbackTable,~] = anova1(feedbackDta, feedbackGrps, 'off');
    unitInfo(u).TrialEpochStats.FeedbackF = [feedbackTable{2,5}, feedbackTable{2,6}];
    
    % Evaluate Unit Firing Between Early and Late Trial Periods
    frTrlPrdsDta = [curUniErlyTrialActivity(inSeqTrialLog & correctTrialLog), curUniLtTrialActivity(inSeqTrialLog & correctTrialLog)]';
    frTrlPrdsGrp = [ones(1,sum(inSeqTrialLog & correctTrialLog)), ones(1,sum(inSeqTrialLog & correctTrialLog))*2]';
    [~,frTrialPrdTable,~] = anova1(frTrlPrdsDta, frTrlPrdsGrp, 'off');
    unitInfo(u).TrialEpochStats.TrialPeriodsF = [frTrialPrdTable{2,5}, frTrialPrdTable{2,6}];
    
    % Evaluate Unit Epoch Activity Modulation by Odor
    [~,frPreTrialOdorTable,~] = anova1(curUniPreTrialActivity(correctTrialLog), trialOdor(correctTrialLog), 'off');
    unitInfo(u).TrialEpochStats.PreTrial.Odor = [frPreTrialOdorTable{2,5}, frPreTrialOdorTable{2,6}];
    [~,frErlyTrialOdorTable,~] = anova1(curUniErlyTrialActivity(correctTrialLog), trialOdor(correctTrialLog), 'off');
    unitInfo(u).TrialEpochStats.EarlyTrial.Odor = [frErlyTrialOdorTable{2,5}, frErlyTrialOdorTable{2,6}];
    [~,frLateTrialOdorTable,~] = anova1(curUniLtTrialActivity(correctTrialLog), trialOdor(correctTrialLog), 'off');
    unitInfo(u).TrialEpochStats.LateTrial.Odor = [frLateTrialOdorTable{2,5}, frLateTrialOdorTable{2,6}];
    [~,frPostTrialOdorTable,~] = anova1(curUniPstTrlActivity(correctTrialLog), trialOdor(correctTrialLog), 'off');
    unitInfo(u).TrialEpochStats.PostTrial.Odor = [frPostTrialOdorTable{2,5}, frPostTrialOdorTable{2,6}];
    [~,frRwdOdorTable,~] = anova1(curUniRwdTrlActivity(correctTrialLog), trialOdor(correctTrialLog), 'off');
    unitInfo(u).TrialEpochStats.Reward.Odor = [frRwdOdorTable{2,5}, frRwdOdorTable{2,6}];
    [~,frErrOdorTable,~] = anova1(curUniErrTrlActivity(~correctTrialLog), trialOdor(~correctTrialLog), 'off');
    unitInfo(u).TrialEpochStats.Error.Odor = [frErrOdorTable{2,5}, frErrOdorTable{2,6}];
    
    % Evaluate Unit Epoch Activity Modulation by Position
    [~,frPreTrialPositionTable,~] = anova1(curUniPreTrialActivity(correctTrialLog), trialPosition(correctTrialLog), 'off');
    unitInfo(u).TrialEpochStats.PreTrial.Position = [frPreTrialPositionTable{2,5}, frPreTrialPositionTable{2,6}];
    [~,frErlyTrialPositionTable,~] = anova1(curUniErlyTrialActivity(correctTrialLog), trialPosition(correctTrialLog), 'off');
    unitInfo(u).TrialEpochStats.EarlyTrial.Position = [frErlyTrialPositionTable{2,5}, frErlyTrialPositionTable{2,6}];
    [~,frLateTrialPositionTable,~] = anova1(curUniLtTrialActivity(correctTrialLog), trialPosition(correctTrialLog), 'off');
    unitInfo(u).TrialEpochStats.LateTrial.Position = [frLateTrialPositionTable{2,5}, frLateTrialPositionTable{2,6}];
    [~,frPostTrialPositionTable,~] = anova1(curUniPstTrlActivity(correctTrialLog), trialPosition(correctTrialLog), 'off');
    unitInfo(u).TrialEpochStats.PostTrial.Position = [frPostTrialPositionTable{2,5}, frPostTrialPositionTable{2,6}];
    [~,frRwdPositionTable,~] = anova1(curUniRwdTrlActivity(correctTrialLog), trialPosition(correctTrialLog), 'off');
    unitInfo(u).TrialEpochStats.Reward.Position = [frRwdPositionTable{2,5}, frRwdPositionTable{2,6}];
    [~,frErrPositionTable,~] = anova1(curUniErrTrlActivity(~correctTrialLog), trialPosition(~correctTrialLog), 'off');
    unitInfo(u).TrialEpochStats.Error.Position = [frErrPositionTable{2,5}, frErrPositionTable{2,6}];
    
    % Evaluate Unit Epoch Activity Modulation by Previous Odor
    [~,frPreTrialPrevOdrTable,~] = anova1(curUniPreTrialActivity(prevOdrLog), prevOdrVect(prevOdrLog), 'off');
    unitInfo(u).TrialEpochStats.PreTrial.PrevOdor = [frPreTrialPrevOdrTable{2,5}, frPreTrialPrevOdrTable{2,6}];
    [~,frErlyTrialPrevOdrTable,~] = anova1(curUniErlyTrialActivity(prevOdrLog), prevOdrVect(prevOdrLog), 'off');
    unitInfo(u).TrialEpochStats.EarlyTrial.PrevOdor = [frErlyTrialPrevOdrTable{2,5}, frErlyTrialPrevOdrTable{2,6}];
    [~,frLateTrialPrevOdrTable,~] = anova1(curUniLtTrialActivity(prevOdrLog), prevOdrVect(prevOdrLog), 'off');
    unitInfo(u).TrialEpochStats.LateTrial.PrevOdor = [frLateTrialPrevOdrTable{2,5}, frLateTrialPrevOdrTable{2,6}];
    [~,frPostTrialPrevOdrTable,~] = anova1(curUniPstTrlActivity(prevOdrLog), prevOdrVect(prevOdrLog), 'off');
    unitInfo(u).TrialEpochStats.PostTrial.PrevOdor = [frPostTrialPrevOdrTable{2,5}, frPostTrialPrevOdrTable{2,6}];
    
    % Evaluate Unit Epoch Activity Modulation by Upcoming Position
    [~,frPreTrialNxtPosTable,~] = anova1(curUniPreTrialActivity(nextPosLog), nextPosVect(nextPosLog), 'off');
    unitInfo(u).TrialEpochStats.PreTrial.NxtPos = [frPreTrialNxtPosTable{2,5}, frPreTrialNxtPosTable{2,6}];
    [~,frErlyTrialNxtPosTable,~] = anova1(curUniErlyTrialActivity(nextPosLog), nextPosVect(nextPosLog), 'off');
    unitInfo(u).TrialEpochStats.EarlyTrial.NxtPos = [frErlyTrialNxtPosTable{2,5}, frErlyTrialNxtPosTable{2,6}];
    [~,frLateTrialNxtPosTable,~] = anova1(curUniLtTrialActivity(nextPosLog), nextPosVect(nextPosLog), 'off');
    unitInfo(u).TrialEpochStats.LateTrial.NxtPos = [frLateTrialNxtPosTable{2,5}, frLateTrialNxtPosTable{2,6}];
    [~,frPostTrialNxtPosTable,~] = anova1(curUniPstTrlActivity(nextPosLog), nextPosVect(nextPosLog), 'off');
    unitInfo(u).TrialEpochStats.PostTrial.NxtPos = [frPostTrialNxtPosTable{2,5}, frPostTrialNxtPosTable{2,6}];  
end   

fprintf('Completed\n');
%% Analysis #2: F-Ratio over time... position vs odor during different periods
fprintf('Starting Analysis 2....');
preEarlyTrialBehavMatrix = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [-0.9 0.6], 'PokeIn');
preEarlyTrialEnsemble = cell2mat(reshape(ExtractTrialData_SM(preEarlyTrialBehavMatrix, ensembleMatrix(:,2:end)), [1,1,length(preEarlyTrialBehavMatrix)])); 

[posFvalsEarlyTrial, ~, posFvalsEarlyTrialZ] = UnitFvalCalcPERM_SM(preEarlyTrialEnsemble(:,:,prevOdrLog), curPosVect(prevOdrLog), slideWindowSize, numPerms);
[odrFvalsEarlyTrial, ~, odrFvalsEarlyTrialZ] = UnitFvalCalcPERM_SM(preEarlyTrialEnsemble(:,:,prevOdrLog), prevOdrVect(prevOdrLog), slideWindowSize, numPerms);

latePostTrialBehavMatrix = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [-0.6 0.6], 'PokeOut');
latePostTrialEnsemble = cell2mat(reshape(ExtractTrialData_SM(latePostTrialBehavMatrix, ensembleMatrix(:,2:end)), [1,1,length(latePostTrialBehavMatrix)])); 

[posFvalsLateTrial, ~, posFvalsLateTrialZ] = UnitFvalCalcPERM_SM(latePostTrialEnsemble(:,:,nextPosLog), nextPosVect(nextPosLog), slideWindowSize, numPerms);
[odrFvalsLateTrial, ~, odrFvalsLateTrialZ] = UnitFvalCalcPERM_SM(latePostTrialEnsemble(:,:,nextPosLog), curOdrVect(nextPosLog), slideWindowSize, numPerms);

for u = 1:length(unitIDs)
    unitInfo(u).InformationContent.EarlyTrial.CurrPosRaw = posFvalsEarlyTrial(:,u);
    unitInfo(u).InformationContent.EarlyTrial.CurrPosZ = posFvalsEarlyTrialZ(:,u);
    unitInfo(u).InformationContent.EarlyTrial.PrevOdorRaw = odrFvalsEarlyTrial(:,u);
    unitInfo(u).InformationContent.EarlyTrial.PrevOdorZ = odrFvalsEarlyTrialZ(:,u);
    
    unitInfo(u).InformationContent.LateTrial.NextPosRaw = posFvalsLateTrial(:,u);
    unitInfo(u).InformationContent.LateTrial.NextPosZ = posFvalsLateTrialZ(:,u);
    unitInfo(u).InformationContent.LateTrial.CurOdorRaw = odrFvalsLateTrial(:,u);
    unitInfo(u).InformationContent.LateTrial.CurOdorZ = odrFvalsLateTrialZ(:,u);
end

fprintf('Completed\n');
%% Analysis #3: Examine Information content by LFP Phase
% This is best done using the epoch extraction script since it give lfp
% phase values.
[earlyUnitEpoch, earlyUnitIDs, earlyLfpEpoch, earlyLfpIDs, ~, ~, earlyTrialInfo] = EpochExtraction_SM('PokeIn', -0.9, 0.6, 'org', 'TiUTr', 'lfpBand', 'All', 'lfpData', 'Phase');

currPos = nan(size(earlyTrialInfo,1),1);
prevOdr = nan(size(earlyTrialInfo,1),1);
for trl = 2:size(earlyTrialInfo,1)
    if earlyTrialInfo(trl,1)==1 && (earlyTrialInfo(trl,3) - earlyTrialInfo(trl-1,3) == 1)
        currPos(trl) = earlyTrialInfo(trl,3);
        prevOdr(trl) = earlyTrialInfo(trl-1,4);
    end
end
earlyEpochEnsmbl = earlyUnitEpoch(:,:,~isnan(currPos));
earlyEpochLFP = earlyLfpEpoch(:,:,~isnan(currPos));
currPos = currPos(~isnan(currPos));
prevOdr = prevOdr(~isnan(prevOdr));     

lfpIDparts = cellfun(@(b)[b(1);b(3)], cellfun(@(a)strsplit(a, '_'), earlyLfpIDs, 'uniformoutput', 0), 'uniformoutput', 0);
lfpIDparts = [lfpIDparts{:}];
bands = unique(lfpIDparts(2,:));
bands(strcmp(bands, 'Raw')) = [];

for uni = 1:length(earlyUnitIDs)
    curUniInfoSpot = strcmp(earlyUnitIDs{uni}, {unitInfo.UnitName});
    curTet = earlyUnitIDs{uni}(1:regexp(earlyUnitIDs{uni}, '-')-1);
    curTetEpochLog = logical(earlyEpochEnsmbl(:,uni,:));
    for band = 1:length(bands)
        curTetLFPepoch = earlyEpochLFP(:,strcmp(bands{band},lfpIDparts(2,:)) & strcmp(curTet,lfpIDparts(1,:)),:);
        curTetSpikePhaseVals = nan(size(curTetEpochLog));
        curTetSpikePhaseVals(curTetEpochLog) = curTetLFPepoch(curTetEpochLog);
        for phase = 1:length(phaseBins)-1
            curTetCurBandPhaseBinSpikes = double(curTetSpikePhaseVals>=phaseBins(phase) & curTetSpikePhaseVals<phaseBins(phase+1));
            [unitInfo(curUniInfoSpot).InformationContentSpikePhase.(bands{band}).(phaseBinLabels{phase}).EarlyTrial.CurrPosRaw, ~,...
                unitInfo(curUniInfoSpot).InformationContentSpikePhase.(bands{band}).(phaseBinLabels{phase}).EarlyTrial.CurrPosZ] = UnitFvalCalcPERM_SM(curTetCurBandPhaseBinSpikes, currPos, slideWindowSize, numPerms);
            
            [unitInfo(curUniInfoSpot).InformationContentSpikePhase.(bands{band}).(phaseBinLabels{phase}).EarlyTrial.PrevOdorRaw, ~,...
                unitInfo(curUniInfoSpot).InformationContentSpikePhase.(bands{band}).(phaseBinLabels{phase}).EarlyTrial.PrevOdorZ] = UnitFvalCalcPERM_SM(curTetCurBandPhaseBinSpikes, prevOdr, slideWindowSize, numPerms);
        end
    end
end          
   
%%%%%%%%% Add in late trial stuff here when there's the luxury of time.
%% Analysis #4: Create Aligned BehavMatrix for the whole trial
fprintf('Starting Analysis 4....');
wholeTrialBehavMatrix = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [-0.9 2], 'PokeIn');
wholeTrialEnsemble = cell2mat(reshape(ExtractTrialData_SM(wholeTrialBehavMatrix, ensembleMatrix(:,2:end)), [1,1,length(wholeTrialBehavMatrix)])); 

eventTimes = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [0 0], 'PokeIn');
trialTimeBins = ExtractTrialData_SM(wholeTrialBehavMatrix, behavMatrix(:,1));
eventTimes = ExtractTrialData_SM(eventTimes, behavMatrix(:,1));
eventTimeBins = cellfun(@(a,b) a-b, trialTimeBins, eventTimes, 'uniformoutput',0);
frstNonMptTrl = find(cellfun(@(a)~isempty(a), eventTimeBins),1, 'first');
eventTimeBins = eventTimeBins{frstNonMptTrl};

for u = 1:length(unitIDs)
    curUniTrlFR = nan(size(wholeTrialEnsemble,1), size(wholeTrialEnsemble,3));
    curUniTrlRstr = cell(1,size(wholeTrialEnsemble,3));
    for trl = 1:size(curUniTrlFR,2)
        curUniTrlFR(:,trl) = conv(wholeTrialEnsemble(:,u,trl), instFRgauss, 'same');
        curUniTrlRstr{trl} = eventTimeBins(logical(wholeTrialEnsemble(:,u,trl)));
    end
    unitInfo(u).WholeTrial.TimeBins = eventTimeBins;
    unitInfo(u).WholeTrial.FiringRate = curUniTrlFR;
    unitInfo(u).WholeTrial.Rasters = curUniTrlRstr;
end
fprintf('Completed\n');
%% Save Analyses
for u = 1:length(unitIDs)
    fprintf('Saving %s UniSum', unitInfo(u).UnitName);
    uniSum = unitInfo(u); %#ok<NASGU>
    save(sprintf('%s_UniSum.mat', unitInfo(u).UnitName), 'uniSum');
    fprintf('... SAVED!\n');
end
