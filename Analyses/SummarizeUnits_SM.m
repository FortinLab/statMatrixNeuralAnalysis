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
numPerms = 10;

%% Analysis #1: Examine Average Evoked activity during trial periods
fprintf('Starting Analysis #1 @%s....', datetime);
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
sansaLog = trialPosition~=1;
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
    
    % Extract Unit Firing Activity SANSA
    curUniPreTrialActivitySANSA = cell2mat(cellfun(@(a)sum(a(:,u)), preTrialEnsemble(sansaLog), 'uniformoutput', 0));       
    curUniErlyTrialActivitySANSA = cell2mat(cellfun(@(a)sum(a(:,u)), earlyTrialEnsemble(sansaLog), 'uniformoutput', 0));  
    curUniLtTrialActivitySANSA = cell2mat(cellfun(@(a)sum(a(:,u)), lateTrialEnsemble(sansaLog), 'uniformoutput', 0));   
    curUniPstTrlActivitySANSA = cell2mat(cellfun(@(a)sum(a(:,u)), postTrialEnsemble(sansaLog), 'uniformoutput', 0));  
    curUniRwdTrlActivitySANSA = cell2mat(cellfun(@(a)sum(a(:,u)), rewardEnsemble(sansaLog), 'uniformoutput', 0));   
    curUniErrTrlActivitySANSA = cell2mat(cellfun(@(a)sum(a(:,u)), errorEnsemble(sansaLog), 'uniformoutput', 0));
        
    unitInfo(u).TrialEpochFRs.PreTrialFRSANSA = [mean(curUniPreTrialActivitySANSA) std(curUniPreTrialActivitySANSA)];
    unitInfo(u).TrialEpochFRs.EarlyTrialFRSANSA = [mean(curUniErlyTrialActivitySANSA) std(curUniErlyTrialActivitySANSA)];
    unitInfo(u).TrialEpochFRs.LateTrialFRSANSA = [mean(curUniLtTrialActivitySANSA) std(curUniLtTrialActivitySANSA)];
    unitInfo(u).TrialEpochFRs.PostTrialFRSANSA = [mean(curUniPstTrlActivitySANSA) std(curUniPstTrlActivitySANSA)];
    unitInfo(u).TrialEpochFRs.RewardFRSANSA = [mean(curUniRwdTrlActivitySANSA(correctTrialLog(sansaLog))) std(curUniRwdTrlActivitySANSA(correctTrialLog(sansaLog)))];
    unitInfo(u).TrialEpochFRs.ErrorFRSANSA = [mean(curUniErrTrlActivitySANSA(~correctTrialLog(sansaLog))) std(curUniErrTrlActivitySANSA(~correctTrialLog(sansaLog)))];
    
    % Evaluate Unit Firing Across Time Periods
    frPeriodDta = [curUniPreTrialActivity(correctTrialLog), curUniErlyTrialActivity(correctTrialLog), curUniLtTrialActivity(correctTrialLog), curUniPstTrlActivity(correctTrialLog)]';
    frPeriodGrp = [ones(1,sum(correctTrialLog)), ones(1,sum(correctTrialLog))*2, ones(1,sum(correctTrialLog))*3, ones(1,sum(correctTrialLog))*4]';
    [~,frPeriodTable,~] = anova1(frPeriodDta, frPeriodGrp, 'off');
    unitInfo(u).TrialEpochStats.TrialEpochsF = [frPeriodTable{2,5}, frPeriodTable{2,6}];
    
    % Evaluate Unit Firing Across Time Periods SANSA
    frPeriodDtaSANSA = [curUniPreTrialActivitySANSA(correctTrialLog(sansaLog)), curUniErlyTrialActivitySANSA(correctTrialLog(sansaLog)), curUniLtTrialActivitySANSA(correctTrialLog(sansaLog)), curUniPstTrlActivitySANSA(correctTrialLog(sansaLog))]';
    frPeriodGrpSANSA = [ones(1,sum(correctTrialLog(sansaLog))), ones(1,sum(correctTrialLog(sansaLog)))*2, ones(1,sum(correctTrialLog(sansaLog)))*3, ones(1,sum(correctTrialLog(sansaLog)))*4]';
    [~,frPeriodTableSANSA,~] = anova1(frPeriodDtaSANSA, frPeriodGrpSANSA, 'off');
    unitInfo(u).TrialEpochStats.TrialEpochsFSANSA = [frPeriodTableSANSA{2,5}, frPeriodTableSANSA{2,6}];
    
    % Evaluate Unit Feeback Properties
    feedbackDta = [curUniRwdTrlActivity(correctTrialLog), curUniErrTrlActivity(~correctTrialLog)];
    feedbackGrps = [ones(1,sum(correctTrialLog)), ones(1,sum(~correctTrialLog))*2];
    [~,feedbackTable,~] = anova1(feedbackDta, feedbackGrps, 'off');
    unitInfo(u).TrialEpochStats.FeedbackF = [feedbackTable{2,5}, feedbackTable{2,6}];
    
    % Evaluate Unit Feeback Properties SANSA
    feedbackDtaSANSA = [curUniRwdTrlActivity(correctTrialLog(sansaLog)), curUniErrTrlActivity(~correctTrialLog(sansaLog))];
    feedbackGrpsSANSA = [ones(1,sum(correctTrialLog(sansaLog))), ones(1,sum(~correctTrialLog(sansaLog)))*2];
    [~,feedbackTableSANSA,~] = anova1(feedbackDtaSANSA, feedbackGrpsSANSA, 'off');
    unitInfo(u).TrialEpochStats.FeedbackFSANSA = [feedbackTableSANSA{2,5}, feedbackTableSANSA{2,6}];
    
    % Evaluate Unit Firing Between Early and Late Trial Periods
    frTrlPrdsDta = [curUniErlyTrialActivity(inSeqTrialLog & correctTrialLog), curUniLtTrialActivity(inSeqTrialLog & correctTrialLog)]';
    frTrlPrdsGrp = [ones(1,sum(inSeqTrialLog & correctTrialLog)), ones(1,sum(inSeqTrialLog & correctTrialLog))*2]';
    [~,frTrialPrdTable,~] = anova1(frTrlPrdsDta, frTrlPrdsGrp, 'off');
    unitInfo(u).TrialEpochStats.TrialPeriodsF = [frTrialPrdTable{2,5}, frTrialPrdTable{2,6}];
    
    % Evaluate Unit Firing Between Early and Late Trial Periods SANSA
    frTrlPrdsDtaSANSA = [curUniErlyTrialActivity(inSeqTrialLog(sansaLog) & correctTrialLog(sansaLog)), curUniLtTrialActivity(inSeqTrialLog(sansaLog) & correctTrialLog(sansaLog))]';
    frTrlPrdsGrpSANSA = [ones(1,sum(inSeqTrialLog(sansaLog) & correctTrialLog(sansaLog))), ones(1,sum(inSeqTrialLog(sansaLog) & correctTrialLog(sansaLog)))*2]';
    [~,frTrialPrdTable,~] = anova1(frTrlPrdsDtaSANSA, frTrlPrdsGrpSANSA, 'off');
    unitInfo(u).TrialEpochStats.TrialPeriodsFSANSA = [frTrialPrdTable{2,5}, frTrialPrdTable{2,6}];
    
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
    
    % Evaluate Unit Epoch Activity Modulation by Odor SANSA
    [~,frPreTrialOdorTableSANSA,~] = anova1(curUniPreTrialActivity(correctTrialLog(sansaLog)), trialOdor(correctTrialLog(sansaLog)), 'off');
    unitInfo(u).TrialEpochStats.PreTrial.OdorSANSA = [frPreTrialOdorTableSANSA{2,5}, frPreTrialOdorTableSANSA{2,6}];
    [~,frErlyTrialOdorTableSANSA,~] = anova1(curUniErlyTrialActivity(correctTrialLog(sansaLog)), trialOdor(correctTrialLog(sansaLog)), 'off');
    unitInfo(u).TrialEpochStats.EarlyTrial.OdorSANSA = [frErlyTrialOdorTableSANSA{2,5}, frErlyTrialOdorTableSANSA{2,6}];
    [~,frLateTrialOdorTableSANSA,~] = anova1(curUniLtTrialActivity(correctTrialLog(sansaLog)), trialOdor(correctTrialLog(sansaLog)), 'off');
    unitInfo(u).TrialEpochStats.LateTrial.OdorSANSA = [frLateTrialOdorTableSANSA{2,5}, frLateTrialOdorTableSANSA{2,6}];
    [~,frPostTrialOdorTableSANSA,~] = anova1(curUniPstTrlActivity(correctTrialLog(sansaLog)), trialOdor(correctTrialLog(sansaLog)), 'off');
    unitInfo(u).TrialEpochStats.PostTrial.OdorSANSA = [frPostTrialOdorTableSANSA{2,5}, frPostTrialOdorTableSANSA{2,6}];
    [~,frRwdOdorTableSANSA,~] = anova1(curUniRwdTrlActivity(correctTrialLog(sansaLog)), trialOdor(correctTrialLog(sansaLog)), 'off');
    unitInfo(u).TrialEpochStats.Reward.OdorSANSA = [frRwdOdorTableSANSA{2,5}, frRwdOdorTableSANSA{2,6}];
    [~,frErrOdorTableSANSA,~] = anova1(curUniErrTrlActivity(~correctTrialLog(sansaLog)), trialOdor(~correctTrialLog(sansaLog)), 'off');
    unitInfo(u).TrialEpochStats.Error.OdorSANSA = [frErrOdorTableSANSA{2,5}, frErrOdorTableSANSA{2,6}];
    
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
    
     % Evaluate Unit Epoch Activity Modulation by Position SANSA
    [~,frPreTrialPositionTableSANSA,~] = anova1(curUniPreTrialActivity(correctTrialLog(sansaLog)), trialPosition(correctTrialLog(sansaLog)), 'off');
    unitInfo(u).TrialEpochStats.PreTrial.PositionSANSA = [frPreTrialPositionTableSANSA{2,5}, frPreTrialPositionTableSANSA{2,6}];
    [~,frErlyTrialPositionTableSANSA,~] = anova1(curUniErlyTrialActivity(correctTrialLog(sansaLog)), trialPosition(correctTrialLog(sansaLog)), 'off');
    unitInfo(u).TrialEpochStats.EarlyTrial.PositionSANSA = [frErlyTrialPositionTableSANSA{2,5}, frErlyTrialPositionTableSANSA{2,6}];
    [~,frLateTrialPositionTableSANSA,~] = anova1(curUniLtTrialActivity(correctTrialLog(sansaLog)), trialPosition(correctTrialLog(sansaLog)), 'off');
    unitInfo(u).TrialEpochStats.LateTrial.PositionSANSA = [frLateTrialPositionTableSANSA{2,5}, frLateTrialPositionTableSANSA{2,6}];
    [~,frPostTrialPositionTableSANSA,~] = anova1(curUniPstTrlActivity(correctTrialLog(sansaLog)), trialPosition(correctTrialLog(sansaLog)), 'off');
    unitInfo(u).TrialEpochStats.PostTrial.PositionSANSA = [frPostTrialPositionTableSANSA{2,5}, frPostTrialPositionTableSANSA{2,6}];
    [~,frRwdPositionTableSANSA,~] = anova1(curUniRwdTrlActivity(correctTrialLog(sansaLog)), trialPosition(correctTrialLog(sansaLog)), 'off');
    unitInfo(u).TrialEpochStats.Reward.PositionSANSA = [frRwdPositionTableSANSA{2,5}, frRwdPositionTableSANSA{2,6}];
    [~,frErrPositionTableSANSA,~] = anova1(curUniErrTrlActivity(~correctTrialLog(sansaLog)), trialPosition(~correctTrialLog(sansaLog)), 'off');
    unitInfo(u).TrialEpochStats.Error.PositionSANSA = [frErrPositionTableSANSA{2,5}, frErrPositionTableSANSA{2,6}];
        
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
    
    % Evaluate Unit Epoch Correlations
    curEpochCorrTable = nan(6);
    curEpochSigTable = nan(6);
    %   Pre-Trial lead
    [peR, peS] = corrcoef(curUniPreTrialActivity(correctTrialLog & inSeqTrialLog), curUniErlyTrialActivity(correctTrialLog & inSeqTrialLog));
    curEpochCorrTable(1,2) = peR(2);
    curEpochSigTable(1,2) = peS(2);
    [plR, plS] = corrcoef(curUniPreTrialActivity(correctTrialLog & inSeqTrialLog), curUniLtTrialActivity(correctTrialLog & inSeqTrialLog));
    curEpochCorrTable(1,3) = plR(2);
    curEpochSigTable(1,3) = plS(2);
    [ppR, ppS] = corrcoef(curUniPreTrialActivity(correctTrialLog & inSeqTrialLog), curUniPstTrlActivity(correctTrialLog & inSeqTrialLog));
    curEpochCorrTable(1,4) = ppR(2);
    curEpochSigTable(1,4) = ppS(2);
    [prR, prS] = corrcoef(curUniPreTrialActivity(correctTrialLog & inSeqTrialLog), curUniRwdTrlActivity(correctTrialLog & inSeqTrialLog));
    curEpochCorrTable(1,5) = prR(2);
    curEpochSigTable(1,5) = prS(2);
    [peR, peS] = corrcoef(curUniPreTrialActivity(~correctTrialLog & inSeqTrialLog), curUniErrTrlActivity(~correctTrialLog & inSeqTrialLog));
    curEpochCorrTable(1,6) = peR(2);
    curEpochSigTable(1,6) = peS(2);
    %   Early-Trial lead
    [elR, elS] = corrcoef(curUniErlyTrialActivity(correctTrialLog & inSeqTrialLog), curUniLtTrialActivity(correctTrialLog & inSeqTrialLog));
    curEpochCorrTable(2,3) = elR(2);
    curEpochSigTable(2,3) = elS(2);
    [epR, epS] = corrcoef(curUniErlyTrialActivity(correctTrialLog & inSeqTrialLog), curUniPstTrlActivity(correctTrialLog & inSeqTrialLog));
    curEpochCorrTable(2,4) = epR(2);
    curEpochSigTable(2,4) = epS(2);
    [erR, erS] = corrcoef(curUniErlyTrialActivity(correctTrialLog & inSeqTrialLog), curUniRwdTrlActivity(correctTrialLog & inSeqTrialLog));
    curEpochCorrTable(2,5) = erR(2);
    curEpochSigTable(2,5) = erS(2);
    [eeR, eeS] = corrcoef(curUniErlyTrialActivity(~correctTrialLog & inSeqTrialLog), curUniErrTrlActivity(~correctTrialLog & inSeqTrialLog));
    curEpochCorrTable(2,6) = eeR(2);
    curEpochSigTable(2,6) = eeS(2);
    %   Late-Trial lead
    [lpR, lpS] = corrcoef(curUniLtTrialActivity(correctTrialLog & inSeqTrialLog), curUniPstTrlActivity(correctTrialLog & inSeqTrialLog));  
    curEpochCorrTable(3,4) = lpR(2);
    curEpochSigTable(3,4) = lpS(2);
    [lrR, lrS] = corrcoef(curUniLtTrialActivity(correctTrialLog & inSeqTrialLog), curUniRwdTrlActivity(correctTrialLog & inSeqTrialLog));
    curEpochCorrTable(3,5) = lrR(2);
    curEpochSigTable(3,5) = lrS(2);
    [leR, leS] = corrcoef(curUniLtTrialActivity(~correctTrialLog & inSeqTrialLog), curUniErrTrlActivity(~correctTrialLog & inSeqTrialLog));
    curEpochCorrTable(3,6) = leR(2);
    curEpochSigTable(3,6) = leS(2);
    %   Post-Trial lead
    [prR, prS] = corrcoef(curUniPstTrlActivity(correctTrialLog & inSeqTrialLog), curUniRwdTrlActivity(correctTrialLog & inSeqTrialLog));
    curEpochCorrTable(4,5) = prR(2);
    curEpochSigTable(4,5) = prS(2);
    [peR, peS] = corrcoef(curUniPstTrlActivity(~correctTrialLog & inSeqTrialLog), curUniErrTrlActivity(~correctTrialLog & inSeqTrialLog));
    curEpochCorrTable(4,6) = peR(2);
    curEpochSigTable(4,6) = peS(2);
    
    unitInfo(u).TrialEpochStats.EpochCorrelations.R = curEpochCorrTable;
    unitInfo(u).TrialEpochStats.EpochCorrelations.P = curEpochSigTable;
%     figure; imagesc(curEpochCorrTable, [-1 1]);
%     set(gca, 'xticklabels', {'Pre-Trial', 'Early-Trial', 'Late-Trial', 'Post-Trial', 'Reward', 'Error'}, 'XTickLabelRotation', 45,...
%         'yticklabels', {'Pre-Trial', 'Early-Trial', 'Late-Trial', 'Post-Trial', 'Reward', 'Error'});
%     hold on;
%     [col,row] = ind2sub([6 6],find(curEpochSigTable<0.0001));
%     for roW = 1:length(row)
%         text(row, col, '*', 'horizontalalignment', 'center', 'fontweight', 'bold', 'fontsize', 50, 'color', 'white')
%     end
%     colormap jet
%     drawnow
end   

fprintf('Completed\n');
%% Analysis #2: F-Ratio over time... position vs odor during different periods
fprintf('Starting Analysis 2 @%s....', datetime);
%%% Poke In Aligned
trialPokeInBehavMatrix = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [-0.9 2.1], 'PokeIn');
trialPokeInEnsemble = cell2mat(reshape(ExtractTrialData_SM(trialPokeInBehavMatrix, ensembleMatrix(:,2:end)), [1,1,length(trialPokeInBehavMatrix)])); 
eventTimes = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [0 0], 'PokeIn');
trialTimeBins = ExtractTrialData_SM(trialPokeInBehavMatrix, behavMatrix(:,1));
eventTimes = ExtractTrialData_SM(eventTimes, behavMatrix(:,1));
eventTimeBins = cellfun(@(a,b) a-b, trialTimeBins, eventTimes, 'uniformoutput',0);
frstNonMptTrl = find(cellfun(@(a)~isempty(a), eventTimeBins),1, 'first');
trlPokeInTimeBins = eventTimeBins{frstNonMptTrl};

% All Correct Trials
[posFvalsPokeInTrial, ~, posFvalsPokeInTrialZ, trialPokeInTimeBins] = UnitFvalCalcPERM_SM(trialPokeInEnsemble(:,:,correctTrialLog), trialPosition(correctTrialLog), slideWindowSize, numPerms, trlPokeInTimeBins);
[odrFvalsPokeInTrial, ~, odrFvalsPokeInTrialZ,~] = UnitFvalCalcPERM_SM(trialPokeInEnsemble(:,:,correctTrialLog), trialOdor(correctTrialLog), slideWindowSize, numPerms, trlPokeInTimeBins);
% Sans Odor A
[posFvalsPokeInSANSATrial, ~, posFvalsPokeInSANSATrialZ, ~] = UnitFvalCalcPERM_SM(trialPokeInEnsemble(:,:,correctTrialLog & sansaLog), trialPosition(correctTrialLog & sansaLog), slideWindowSize, numPerms, trlPokeInTimeBins);
[odrFvalsPokeInSANSATrial, ~, odrFvalsPokeInSANSATrialZ, ~] = UnitFvalCalcPERM_SM(trialPokeInEnsemble(:,:,correctTrialLog & sansaLog), trialOdor(correctTrialLog & sansaLog), slideWindowSize, numPerms, trlPokeInTimeBins);
[prevOdrFvalsPokeIn, ~, prevOdrFvalsPokeInZ, ~] = UnitFvalCalcPERM_SM(trialPokeInEnsemble(:,:,prevOdrLog), prevOdrVect(prevOdrLog), slideWindowSize, numPerms, trlPokeInTimeBins);
% Sans Odor A All In Sequence
[posFvalsPokeInSANSAaisTrial, ~, posFvalsPokeInSANSAaisTrialZ, ~] = UnitFvalCalcPERM_SM(trialPokeInEnsemble(:,:,correctTrialLog & sansaLog & inSeqTrialLog), trialPosition(correctTrialLog & sansaLog & inSeqTrialLog), slideWindowSize, numPerms, trlPokeInTimeBins);

%%% Poke Out Aligned
trialPokeOutBehavMatrix = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [-1.9 1.1], 'PokeOut');
trialPokeOutEnsemble = cell2mat(reshape(ExtractTrialData_SM(trialPokeOutBehavMatrix, ensembleMatrix(:,2:end)), [1,1,length(trialPokeOutBehavMatrix)])); 
eventTimes = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [0 0], 'PokeOut');
trialTimeBins = ExtractTrialData_SM(trialPokeOutBehavMatrix, behavMatrix(:,1));
eventTimes = ExtractTrialData_SM(eventTimes, behavMatrix(:,1));
eventTimeBins = cellfun(@(a,b) a-b, trialTimeBins, eventTimes, 'uniformoutput',0);
frstNonMptTrl = find(cellfun(@(a)~isempty(a), eventTimeBins),1, 'first');
trlPokeOutTimeBins = eventTimeBins{frstNonMptTrl};

% All Correct Trials
[posFvalsPokeOutTrial, ~, posFvalsPokeOutTrialZ, trialPokeOutTimeBins] = UnitFvalCalcPERM_SM(trialPokeOutEnsemble(:,:,correctTrialLog), trialPosition(correctTrialLog), slideWindowSize, numPerms, trlPokeOutTimeBins);
[odrFvalsPokeOutTrial, ~, odrFvalsPokeOutTrialZ,~] = UnitFvalCalcPERM_SM(trialPokeOutEnsemble(:,:,correctTrialLog), trialOdor(correctTrialLog), slideWindowSize, numPerms, trlPokeOutTimeBins);
% Sans Odor A
[posFvalsPokeOutSANSATrial, ~, posFvalsPokeOutSANSATrialZ, ~] = UnitFvalCalcPERM_SM(trialPokeOutEnsemble(:,:,correctTrialLog & sansaLog), trialPosition(correctTrialLog & sansaLog), slideWindowSize, numPerms, trlPokeOutTimeBins);
[odrFvalsPokeOutSANSATrial, ~, odrFvalsPokeOutSANSATrialZ,~] = UnitFvalCalcPERM_SM(trialPokeOutEnsemble(:,:,correctTrialLog & sansaLog), trialOdor(correctTrialLog & sansaLog), slideWindowSize, numPerms, trlPokeOutTimeBins);
[prevOdrFvalsPokeOut, ~, prevOdrFvalsPokeOutZ, ~] = UnitFvalCalcPERM_SM(trialPokeOutEnsemble(:,:,prevOdrLog), prevOdrVect(prevOdrLog), slideWindowSize, numPerms, trlPokeOutTimeBins);
% Sans Odor A All In Sequence
[posFvalsPokeOutSANSAaisTrial, ~, posFvalsPokeOutSANSAaisTrialZ, ~] = UnitFvalCalcPERM_SM(trialPokeOutEnsemble(:,:,correctTrialLog & sansaLog & inSeqTrialLog), trialPosition(correctTrialLog & sansaLog & inSeqTrialLog), slideWindowSize, numPerms, trlPokeOutTimeBins);

%%% Compile Everything into UnitInfo structures
for u = 1:length(unitIDs)
    % Poke In Aligned
    % All Correct Trials
    unitInfo(u).InformationContent.TrialPokeIn.TimeBins = trialPokeInTimeBins;
    unitInfo(u).InformationContent.TrialPokeIn.PosRaw = posFvalsPokeInTrial(:,u);
    unitInfo(u).InformationContent.TrialPokeIn.PosZ = posFvalsPokeInTrialZ(:,u);
    unitInfo(u).InformationContent.TrialPokeIn.OdorRaw = odrFvalsPokeInTrial(:,u);
    unitInfo(u).InformationContent.TrialPokeIn.OdorZ = odrFvalsPokeInTrialZ(:,u);
    % All correct trials Sans Odor A
    unitInfo(u).InformationContent.TrialPokeInSANSA.TimeBins = trialPokeInTimeBins;
    unitInfo(u).InformationContent.TrialPokeInSANSA.PosRaw = posFvalsPokeInSANSATrial(:,u);
    unitInfo(u).InformationContent.TrialPokeInSANSA.PosZ = posFvalsPokeInSANSATrialZ(:,u);
    unitInfo(u).InformationContent.TrialPokeInSANSA.OdorRaw = odrFvalsPokeInSANSATrial(:,u);
    unitInfo(u).InformationContent.TrialPokeInSANSA.OdorZ = odrFvalsPokeInSANSATrialZ(:,u);  
    unitInfo(u).InformationContent.TrialPokeInSANSA.PrevOdorRaw = prevOdrFvalsPokeIn(:,u);
    unitInfo(u).InformationContent.TrialPokeInSANSA.PrevOdorZ = prevOdrFvalsPokeInZ(:,u);
    % All In Sequence Sans Odor A
    unitInfo(u).InformationContent.TrialPokeInSANSAais.TimeBins = trialPokeInTimeBins;
    unitInfo(u).InformationContent.TrialPokeInSANSAais.PosRaw = posFvalsPokeInSANSAaisTrial(:,u);
    unitInfo(u).InformationContent.TrialPokeInSANSAais.PosZ = posFvalsPokeInSANSAaisTrialZ(:,u);
    
    % Poke Out Aligned
    % All Correct Trials
    unitInfo(u).InformationContent.TrialPokeOut.TimeBins = trialPokeOutTimeBins;
    unitInfo(u).InformationContent.TrialPokeOut.PosRaw = posFvalsPokeOutTrial(:,u);
    unitInfo(u).InformationContent.TrialPokeOut.PosZ = posFvalsPokeOutTrialZ(:,u);
    unitInfo(u).InformationContent.TrialPokeOut.OdorRaw = odrFvalsPokeOutTrial(:,u);
    unitInfo(u).InformationContent.TrialPokeOut.OdorZ = odrFvalsPokeOutTrialZ(:,u);
    % All correct trials Sans Odor A
    unitInfo(u).InformationContent.TrialPokeOutSANSA.TimeBins = trialPokeOutTimeBins;
    unitInfo(u).InformationContent.TrialPokeOutSANSA.PosRaw = posFvalsPokeOutSANSATrial(:,u);
    unitInfo(u).InformationContent.TrialPokeOutSANSA.PosZ = posFvalsPokeOutSANSATrialZ(:,u);
    unitInfo(u).InformationContent.TrialPokeOutSANSA.OdorRaw = odrFvalsPokeOutSANSATrial(:,u);
    unitInfo(u).InformationContent.TrialPokeOutSANSA.OdorZ = odrFvalsPokeOutSANSATrialZ(:,u);
    unitInfo(u).InformationContent.TrialPokeOutSANSA.PrevOdorRaw = prevOdrFvalsPokeOut(:,u);
    unitInfo(u).InformationContent.TrialPokeOutSANSA.PrevOdorZ = prevOdrFvalsPokeOutZ(:,u);
    % All In Sequence Sans Odor A
    unitInfo(u).InformationContent.TrialPokeOutSANSAais.TimeBins = trialPokeOutTimeBins;
    unitInfo(u).InformationContent.TrialPokeOutSANSAais.PosRaw = posFvalsPokeOutSANSAaisTrial(:,u);
    unitInfo(u).InformationContent.TrialPokeOutSANSAais.PosZ = posFvalsPokeOutSANSAaisTrialZ(:,u);
end

fprintf('Completed\n');
%% Analysis #3: Examine Information content by LFP Phase
fprintf('Starting Analysis 3 @%s....', datetime);
% This is best done using the epoch extraction script since it give lfp
% phase values.
[earlyUnitEpoch, earlyUnitIDs, earlyLfpEpoch, earlyLfpIDs, ~, rlyTimeBins, earlyTrialInfo] = EpochExtraction_SM('PokeIn', -0.9, 2.1, 'org', 'TiUTr', 'lfpBand', 'All', 'lfpData', 'Phase');

currPos = nan(size(earlyTrialInfo,1),1);
currOdr = nan(size(earlyTrialInfo,1),1);
prevOdr = nan(size(earlyTrialInfo,1),1);
for trl = 2:size(earlyTrialInfo,1)
    if earlyTrialInfo(trl,1)==1 && (earlyTrialInfo(trl,3) - earlyTrialInfo(trl-1,3) == 1)
        currPos(trl) = earlyTrialInfo(trl,3);
        currOdr(trl) = earlyTrialInfo(trl,4);
        prevOdr(trl) = earlyTrialInfo(trl-1,4);
    end
end
earlyEpochEnsmbl = earlyUnitEpoch(:,:,~isnan(currPos));
earlyEpochLFP = earlyLfpEpoch(:,:,~isnan(currPos));
currPos = currPos(~isnan(currPos));
currOdr = currOdr(~isnan(currOdr));
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
                unitInfo(curUniInfoSpot).InformationContentSpikePhase.(bands{band}).(phaseBinLabels{phase}).EarlyTrial.CurrPosZ,earlyTimeBins] = UnitFvalCalcPERM_SM(curTetCurBandPhaseBinSpikes, currPos, slideWindowSize, numPerms, rlyTimeBins);
            
            [unitInfo(curUniInfoSpot).InformationContentSpikePhase.(bands{band}).(phaseBinLabels{phase}).EarlyTrial.CurrOdorRaw, ~,...
                unitInfo(curUniInfoSpot).InformationContentSpikePhase.(bands{band}).(phaseBinLabels{phase}).EarlyTrial.CurrOdorZ,~] = UnitFvalCalcPERM_SM(curTetCurBandPhaseBinSpikes, currOdr, slideWindowSize, numPerms, rlyTimeBins);
            
            [unitInfo(curUniInfoSpot).InformationContentSpikePhase.(bands{band}).(phaseBinLabels{phase}).EarlyTrial.PrevOdorRaw, ~,...
                unitInfo(curUniInfoSpot).InformationContentSpikePhase.(bands{band}).(phaseBinLabels{phase}).EarlyTrial.PrevOdorZ,~] = UnitFvalCalcPERM_SM(curTetCurBandPhaseBinSpikes, prevOdr, slideWindowSize, numPerms, rlyTimeBins);
            unitInfo(curUniInfoSpot).InformationContentSpikePhase.(bands{band}).(phaseBinLabels{phase}).EarlyTrial.TimeBins = earlyTimeBins;
        end
    end
end          
   
%%%%%%%% Add in late trial stuff here when there's the luxury of time.
fprintf('Completed\n');
%% Analysis #4: Create Aligned BehavMatrix for the whole trial
fprintf('Starting Analysis 4 @%s....', datetime);
wholeTrialPokeInBehavMatrix = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [-0.8 2], 'PokeIn');
wholeTrialPokeInEnsemble = cell2mat(reshape(ExtractTrialData_SM(wholeTrialPokeInBehavMatrix, ensembleMatrix(:,2:end)), [1,1,length(wholeTrialPokeInBehavMatrix)])); 

eventTimes = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [0 0], 'PokeIn');
trialTimeBins = ExtractTrialData_SM(wholeTrialPokeInBehavMatrix, behavMatrix(:,1));
eventTimes = ExtractTrialData_SM(eventTimes, behavMatrix(:,1));
eventTimeBins = cellfun(@(a,b) a-b, trialTimeBins, eventTimes, 'uniformoutput',0);
frstNonMptTrl = find(cellfun(@(a)~isempty(a), eventTimeBins),1, 'first');
eventTimeBins = eventTimeBins{frstNonMptTrl};

for u = 1:length(unitIDs)
    curUniTrlFR = nan(size(wholeTrialPokeInEnsemble,1), size(wholeTrialPokeInEnsemble,3));
    curUniTrlRstr = cell(1,size(wholeTrialPokeInEnsemble,3));
    for trl = 1:size(curUniTrlFR,2)
        curUniTrlFR(:,trl) = conv(wholeTrialPokeInEnsemble(:,u,trl), instFRgauss, 'same');
        curUniTrlRstr{trl} = eventTimeBins(logical(wholeTrialPokeInEnsemble(:,u,trl)));
    end
    unitInfo(u).WholeTrial.PokeIn.TimeBins = eventTimeBins;
    unitInfo(u).WholeTrial.PokeIn.FiringRate = curUniTrlFR;
    unitInfo(u).WholeTrial.PokeIn.Rasters = curUniTrlRstr;
end

wholeTrialPokeOutBehavMatrix = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [-1.8 1], 'PokeOut');
wholeTrialPokeOutEnsemble = cell2mat(reshape(ExtractTrialData_SM(wholeTrialPokeOutBehavMatrix, ensembleMatrix(:,2:end)), [1,1,length(wholeTrialPokeOutBehavMatrix)])); 

eventTimes = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [0 0], 'PokeOut');
trialTimeBins = ExtractTrialData_SM(wholeTrialPokeOutBehavMatrix, behavMatrix(:,1));
eventTimes = ExtractTrialData_SM(eventTimes, behavMatrix(:,1));
eventTimeBins = cellfun(@(a,b) a-b, trialTimeBins, eventTimes, 'uniformoutput',0);
frstNonMptTrl = find(cellfun(@(a)~isempty(a), eventTimeBins),1, 'first');
eventTimeBins = eventTimeBins{frstNonMptTrl};

for u = 1:length(unitIDs)
    curUniTrlFR = nan(size(wholeTrialPokeOutEnsemble,1), size(wholeTrialPokeOutEnsemble,3));
    curUniTrlRstr = cell(1,size(wholeTrialPokeOutEnsemble,3));
    for trl = 1:size(curUniTrlFR,2)
        curUniTrlFR(:,trl) = conv(wholeTrialPokeOutEnsemble(:,u,trl), instFRgauss, 'same');
        curUniTrlRstr{trl} = eventTimeBins(logical(wholeTrialPokeOutEnsemble(:,u,trl)));
    end
    unitInfo(u).WholeTrial.PokeOut.TimeBins = eventTimeBins;
    unitInfo(u).WholeTrial.PokeOut.FiringRate = curUniTrlFR;
    unitInfo(u).WholeTrial.PokeOut.Rasters = curUniTrlRstr;
end
fprintf('Completed\n');
%% Save Analyses
dirParts = strsplit(cd, '\');
for u = 1:length(unitIDs)
    fprintf('Saving %s UniSum', unitInfo(u).UnitName);
    uniSum = unitInfo(u); %#ok<NASGU>
    save(sprintf('%s_%s_UniSum.mat', dirParts{end}, unitInfo(u).UnitName), 'uniSum');
    fprintf('... SAVED!\n');
end
PlotUniSum_SM
