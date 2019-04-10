clear all
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
uniSumLog = false(size(matFiles));
for fl = 1:length(matFiles)
    variableInfo = who('-file', matFiles{fl});
    if sum(ismember(variableInfo, 'uniSum'))==1
        uniSumLog(fl) = true;
    end
end
if sum(uniSumLog) == 0
    error('No Unit Summary files located, check directory');
end

uniSumFiles = matFiles(uniSumLog);
for u = 1:length(uniSumFiles)
    load(uniSumFiles{u});
    if ~isfield(uniSum, 'InformationContentSpikePhase')     %******** REMEMBER TO REMOVE ONCE THIS FIELD IS THERE FOR ALL DATA
        uniSum.InformationContentSpikePhase = [];
    end
    uniSumData(u) = uniSum;
end

%% Select tetrodes
tetIDs = cellfun(@(a)a(1),cellfun(@(a)strsplit(a, '-'), {uniSumData.UnitName}, 'uniformoutput', 0));
postLog = strcmp(tetIDs, 'T19') | strcmp(tetIDs, 'T18') | strcmp(tetIDs, 'T17') | strcmp(tetIDs, 'T20') | strcmp(tetIDs, 'T21');
midLog = strcmp(tetIDs, 'T13') | strcmp(tetIDs, 'T12') | strcmp(tetIDs, 'T24') | strcmp(tetIDs, 'T2');
antLog = strcmp(tetIDs, 'T10') | strcmp(tetIDs, 'T9') | strcmp(tetIDs, 'T6') | strcmp(tetIDs, 'T5');

% uniSumData(~postLog) = [];
% uniSumData(~midLog) = [];
% uniSumData(~antLog) = [];

%% Analysis #1: Evaluate population for firing rate differences
newCritF = 0.05/23;
newCritR = 0.05/14;
trialEpochStats = [uniSumData.TrialEpochStats];

trialEpochF = cell2mat(cellfun(@(a){a'},{trialEpochStats.TrialEpochsFSANSA}));
feedbackF = cell2mat(cellfun(@(a){a'},{trialEpochStats.FeedbackFSANSA}));
trialPeriodF = cell2mat(cellfun(@(a){a'},{trialEpochStats.TrialPeriodsFSANSA}));

barVals = [mean(trialEpochF(2,:)<newCritF)*100;...                              % Percentage of neurons with significant trial epoch modulation (pre- vs early- vs late- vs post- vs reward- vs error- periods)
    mean(trialPeriodF(2,:)<newCritF)*100;...                                     % Percentage of neurons significnatly distinguishing the early and late trial periods
    mean(feedbackF(2,:)<newCritF)*100];                                         % Percentage of neurons significantly distinguishing Reward and Error periods
figure;
subplot(3,2,1)
bar(barVals);
set(gca, 'xtick', 1:3, 'xticklabels', {'All Epochs', 'Early vs Late Trial', 'Reward vs Error'}, 'xticklabelrotation', 45, 'ylim', [0 100]);
title('Overall');
ylabel('Percentage of neurons');

trialEpochFRs = [uniSumData.TrialEpochFRs];
uniFRs = [cell2mat(cellfun(@(a){a(1)},{trialEpochFRs.PreTrialFRSANSA}));...
    cell2mat(cellfun(@(a){a(1)},{trialEpochFRs.EarlyTrialFRSANSA}));...
    cell2mat(cellfun(@(a){a(1)},{trialEpochFRs.LateTrialFRSANSA}));...
    cell2mat(cellfun(@(a){a(1)},{trialEpochFRs.PostTrialFRSANSA}));...
    cell2mat(cellfun(@(a){a(1)},{trialEpochFRs.RewardFRSANSA}));...
    cell2mat(cellfun(@(a){a(1)},{trialEpochFRs.ErrorFRSANSA}))];
maxFR = max(uniFRs);
maxFRepoch = nan(size(maxFR));
for u = 1:length(maxFR)
    maxFRepoch(u) = find(uniFRs(:,u)==maxFR(u));
end
subplot(3,2,2)
pie([mean(~(trialEpochF(2,:)<newCritF))*100,...
    mean(maxFRepoch==1 & trialEpochF(2,:)<newCritF)*100,...
    mean(maxFRepoch==2 & trialEpochF(2,:)<newCritF)*100,...
    mean(maxFRepoch==3 & trialEpochF(2,:)<newCritF)*100,...
    mean(maxFRepoch==4 & trialEpochF(2,:)<newCritF)*100,...
    mean(maxFRepoch==5 & trialEpochF(2,:)<newCritF)*100,...
    mean(maxFRepoch==6 & trialEpochF(2,:)<newCritF)*100],...
    true(1,7),...
    {'NR', 'Pre', 'Early', 'Late', 'Post', 'Rwd', 'Err'});


% Select ONLY those modulated by trial epoch
allUniSumData = uniSumData;
uniSumData(~(trialEpochF(2,:)<newCritF)) = [];
maxFRepoch(~(trialEpochF(2,:)<newCritF)) = [];

% subplot(3,2,2)
% barVals = [mean(maxFRepoch==1)*100,...
%     mean(maxFRepoch==2)*100,...
%     mean(maxFRepoch==3)*100,...
%     mean(maxFRepoch==4)*100,...
%     mean(maxFRepoch==5)*100,...
%     mean(maxFRepoch==6)*100];
% bar(barVals);
% set(gca, 'xtick', 1:6, 'xticklabels', {'Pre', 'Early', 'Late', 'Post', 'Rwd', 'Err'}, 'ylim', [0 100], 'xticklabelrotation', 45, 'xlim', [0 7]);
% title('Max Firing Rate Epochs');

trialPeriods = {'PreTrial', 'EarlyTrial', 'LateTrial', 'PostTrial','Reward','Error'};
trlModMtxPos = nan(length(trialPeriods));
trlModMtxOdr = nan(length(trialPeriods));
trlModEpoch = nan(2,length(trialPeriods));
for prd1 = 1:length(trialPeriods)
    curPeriod = [trialEpochStats.(trialPeriods{prd1})];
    curPrdPos = cell2mat(cellfun(@(a){a'},{curPeriod.PositionSANSA}));
    trlModEpoch(1,prd1) = mean(curPrdPos(2,:)<newCritF)*100;
    curPrdOdor = cell2mat(cellfun(@(a){a'},{curPeriod.OdorSANSA}));
    trlModEpoch(2,prd1) = mean(curPrdOdor(2,:)<newCritF)*100;
    for prd2 = 1:length(trialPeriods)
        curPrdLog = maxFRepoch==prd2;
        trlModMtxOdr(prd2,prd1) = mean(curPrdOdor(2,curPrdLog)<newCritF)*100;
        trlModMtxPos(prd2,prd1) = mean(curPrdPos(2,curPrdLog)<newCritF)*100;
    end
end

subplot(3,2,3)
imagesc(1:length(trialPeriods), 1:length(trialPeriods),trlModMtxPos, [0 65]);
set(gca, 'xtick', 1:length(trialPeriods), 'xticklabel', trialPeriods, 'xticklabelrotation', 45, 'ytick', 1:length(trialPeriods), 'yticklabels', trialPeriods, 'dataaspectratio', [1 1 1]);
xlabel('Modulation Period');
ylabel('Max FR Period');
title('Position Modulation');

subplot(3,2,4)
imagesc(1:length(trialPeriods), 1:length(trialPeriods),trlModMtxOdr, [0 65]);
set(gca, 'xtick', 1:length(trialPeriods), 'xticklabel', trialPeriods, 'xticklabelrotation', 45, 'ytick', 1:length(trialPeriods), 'yticklabels', trialPeriods, 'dataaspectratio', [1 1 1]);
xlabel('Modulation Period');
ylabel('Max FR Period');
title('Odor Modulation');

posPrcntBars = subplot(3,2,5)
bar(1:length(trialPeriods),trlModEpoch(1,:));
set(gca, 'xlim', [0 length(trialPeriods)+1], 'xtick', 1:length(trialPeriods), 'xticklabels', trialPeriods, 'xticklabelrotation', 45);
title('Position');
ylabel('Percent Cells Modulated');

odrPrcntBars = subplot(3,2,6);
bar(1:length(trialPeriods), trlModEpoch(2,:));
set(gca, 'xlim', [0 length(trialPeriods)+1], 'xtick', 1:length(trialPeriods), 'xticklabels', trialPeriods, 'xticklabelrotation', 45);
title('Odor');
ylabel('Percent Cells Modulated');

linkaxes([posPrcntBars, odrPrcntBars], 'y');

%% Analysis #2: Create Ensemble Normalized Firing Rate Plot
pokeDurs = uniSumData(1).TrialInfo.PokeDuration;
inSeqPokeDurs = pokeDurs(logical(uniSumData(1).TrialInfo.InSeqTrialLog) & logical(uniSumData(1).TrialInfo.CorrectTrialLog));
meanPokeDur = mean(inSeqPokeDurs);
wholeTrialData = [uniSumData.WholeTrial];

wholeTrialDataPI = [wholeTrialData.PokeIn];
normUniFRsPI = nan(length(wholeTrialDataPI), length(wholeTrialDataPI(1).TimeBins));
peakFRsPI = nan(length(wholeTrialDataPI),1);
for u = 1:length(wholeTrialDataPI)
    curUniFRTrls = wholeTrialDataPI(u).FiringRate;
    curUniFR = mean(curUniFRTrls(:,uniSumData(u).TrialInfo.Position~=1 & uniSumData(u).TrialInfo.CorrectTrialLog & uniSumData(u).TrialInfo.InSeqTrialLog),2);
    normUniFRsPI(u,:) = curUniFR/max(curUniFR);
    peakFRsPI(u) = find(curUniFR==max(curUniFR),1,'first');
end
erlyPosSortingMtxPI = [peakFRsPI, normUniFRsPI];
erlyPosSortedMtxPI = sortrows(erlyPosSortingMtxPI);
sortedFRsPI = erlyPosSortedMtxPI(:,2:end);

piTimeBins = wholeTrialDataPI(1).TimeBins;
% validPItimes = piTimeBins<0.5;
validPItimes = piTimeBins<1000;
piTimeBins(~validPItimes) = [];

figure;
subplot(1,2,1)
imagesc(piTimeBins,1:length(wholeTrialDataPI),sortedFRsPI(:,validPItimes), [0 1]);
title('PFC neurons are active throughout the trial period');
% set(gca, 'xtick', [0 meanPokeDur], 'xticklabel', {'PokeIn', 'Mean PokeOut'}, 'yticklabel', []);
set(gca, 'yticklabel', []);
ylabel([{'PFC Neurons'}; {'Organized by Peak Firing Rate'}]);
xlabel('Time Relative to Trial Start (s)');
hold on;
line([0 0], get(gca, 'ylim'), 'linestyle', '--', 'color', 'white', 'linewidth', 2);
line([meanPokeDur meanPokeDur], get(gca, 'ylim'), 'linestyle', '--', 'color', 'white', 'linewidth', 2);
    

wholeTrialDataPO = [wholeTrialData.PokeOut];
normUniFRsPO = nan(length(wholeTrialDataPO), length(wholeTrialDataPO(1).TimeBins));
peakFRsPO = nan(length(wholeTrialDataPO),1);
for u = 1:length(wholeTrialDataPO)
    curUniFR = mean(wholeTrialDataPO(u).FiringRate,2);
    normUniFRsPO(u,:) = curUniFR/max(curUniFR);
    peakFRsPO(u) = find(curUniFR==max(curUniFR),1,'first');
end
erlyPosSortingMtxPO = [peakFRsPO, normUniFRsPO];
erlyPosSortedMtxPO = sortrows(erlyPosSortingMtxPO);
sortedFRsPO = erlyPosSortedMtxPO(:,2:end);

poTimeBins = wholeTrialDataPO(1).TimeBins;
% validPOtimes = poTimeBins>-0.5;
validPOtimes = poTimeBins>-1000;
poTimeBins(~validPOtimes) = [];

subplot(1,2,2)
imagesc(poTimeBins,1:length(wholeTrialDataPI),sortedFRsPO(:,validPOtimes), [0 1]);
title('PFC neurons are active throughout the trial period');
% set(gca, 'xtick', [meanPokeDur*-1 0], 'xticklabel', {'Mean PokeIn', 'PokeOut'}, 'yticklabel', []);
set(gca, 'yticklabel', []);
ylabel([{'PFC Neurons'}; {'Organized by Peak Firing Rate'}]);
xlabel('Time Relative to Trial End (s)');
hold on;
line([0 0], get(gca, 'ylim'), 'linestyle', '--', 'color', 'white', 'linewidth', 2);
line([meanPokeDur*-1 meanPokeDur*-1], get(gca, 'ylim'), 'linestyle', '--', 'color', 'white', 'linewidth', 2);
    
    
%% Analysis #3: Create Ensemble Normalized Information Plot
ic = [uniSumData.InformationContent];
pokeInIC = [ic.TrialPokeInSANSA];
pokeOutIC = [ic.TrialPokeOutSANSA];

% Early Trial
figure;
pokeInPosZ = [pokeInIC.PosZ];
pokeInOdrZ = [pokeInIC.OdorZ];
pokeInPrevOdrZ = [pokeInIC.PrevOdorZ];
normErlTrlPosZ = nan(size(pokeInPosZ));
normErlTrlOdrZ = nan(size(pokeInOdrZ));
normErlTrlPrevOdrZ = nan(size(pokeInPrevOdrZ));
peakPosZ = nan(1,size(pokeInPosZ,2));
peakOdrZ = nan(1,size(pokeInOdrZ,2));
peakPrevOdrZ = nan(1,size(pokeInPrevOdrZ,2));
for u = 1:size(pokeInPosZ,2);
    curUniPosIC = pokeInPosZ(:,u);
    maxPosIC = max(curUniPosIC);
    peakPosZ(u) = find(curUniPosIC==maxPosIC,1,'first');
    normErlTrlPosZ(:,u) = curUniPosIC./maxPosIC;
    
    curUniOdrIC = pokeInOdrZ(:,u);
    maxOdrIC = max(curUniOdrIC);
    peakOdrZ(u) = find(curUniOdrIC==maxOdrIC,1,'first');
    normErlTrlOdrZ(:,u) = curUniOdrIC./maxOdrIC;
    
    curUniPrevOdrIC = pokeInPrevOdrZ(:,u);
    maxPrevOdrIC = max(curUniPrevOdrIC);
    peakPrevOdrZ(u) = find(curUniPrevOdrIC==maxPrevOdrIC,1,'first');
    normErlTrlPrevOdrZ(:,u) = curUniPrevOdrIC./maxPrevOdrIC;
end
% Position Sorted by Position
erlyPosSortingMtx = [peakPosZ;1:length(peakPosZ);normErlTrlPosZ]';
erlyPosSortedMtx = sortrows(erlyPosSortingMtx);
erlyPosSortedIC = erlyPosSortedMtx(:,3:end);
subplot(3,3,1);
imagesc(pokeInIC(1).TimeBins,1:length(wholeTrialData),erlyPosSortedIC, [0 1]);
title('Position Information');
set(gca, 'xtick', 0, 'xticklabel', 'PokeIn', 'yticklabel', []);
ylabel([{'PFC Neurons'}; {'Organized by Peak F-Value'}]);
xlabel('Time Relative to Trial Start (s)');
hold on;
line([0 0], get(gca, 'ylim'), 'linestyle', '--', 'color', 'white', 'linewidth', 2);
% Position Sorted by Current Odor
erlyPosODRSortingMtx = [peakOdrZ;1:length(peakPosZ);normErlTrlPosZ]';
erlyPosODRSortedMtx = sortrows(erlyPosODRSortingMtx);
erlyPosODRSortedIC = erlyPosODRSortedMtx(:,3:end);
subplot(3,3,2);
imagesc(pokeInIC(1).TimeBins,1:length(wholeTrialData),erlyPosODRSortedIC, [0 1]);
title('Position Sorted by Current Odor');
set(gca, 'xtick', 0, 'xticklabel', 'PokeIn', 'yticklabel', []);
ylabel([{'PFC Neurons'}; {'Organized by Peak F-Value'}]);
xlabel('Time Relative to Trial Start (s)');
hold on;
line([0 0], get(gca, 'ylim'), 'linestyle', '--', 'color', 'white', 'linewidth', 2);
% Position Sorted by Previous Odor
erlyPosPREVODRSortingMtx = [peakPrevOdrZ;1:length(peakPrevOdrZ);normErlTrlPosZ]';
erlyPosPREVODRSortedMtx = sortrows(erlyPosPREVODRSortingMtx);
erlyPosPREVODRSortedIC = erlyPosPREVODRSortedMtx(:,3:end);
subplot(3,3,3);
imagesc(pokeInIC(1).TimeBins,1:length(wholeTrialData),erlyPosPREVODRSortedIC, [0 1]);
title('Position Sorted by Previous Odor');
set(gca, 'xtick', 0, 'xticklabel', 'PokeIn', 'yticklabel', []);
ylabel([{'PFC Neurons'}; {'Organized by Peak F-Value'}]);
xlabel('Time Relative to Trial Start (s)');
hold on;
line([0 0], get(gca, 'ylim'), 'linestyle', '--', 'color', 'white', 'linewidth', 2);

% Current Odor Sorted by Position
erlyOdrPOSSortingMtx = [peakPosZ;1:length(peakPosZ);normErlTrlOdrZ]';
erlyOdrPOSSortedMtx = sortrows(erlyOdrPOSSortingMtx);
erlyOdrPOSSortedIC = erlyOdrPOSSortedMtx(:,3:end);
subplot(3,3,4);
imagesc(pokeInIC(1).TimeBins,1:length(wholeTrialData),erlyOdrPOSSortedIC, [0 1]);
title('Current Odor Sorted by Position');
set(gca, 'xtick', 0, 'xticklabel', 'PokeIn', 'yticklabel', []);
ylabel([{'PFC Neurons'}; {'Organized by Peak F-Value'}]);
xlabel('Time Relative to Trial Start (s)');
hold on;
line([0 0], get(gca, 'ylim'), 'linestyle', '--', 'color', 'white', 'linewidth', 2);
% Odor Sorted by Odor
erlyOdrSortingMtx = [peakOdrZ;1:length(peakOdrZ);normErlTrlOdrZ]';
erlyOdrSortedMtx = sortrows(erlyOdrSortingMtx);
erlyOdrSortedIC = erlyOdrSortedMtx(:,3:end);
subplot(3,3,5);
imagesc(pokeInIC(1).TimeBins,1:length(wholeTrialData),erlyOdrSortedIC, [0 1]);
title('Current Odor Information');
set(gca, 'xtick', 0, 'xticklabel', 'PokeIn', 'yticklabel', []);
ylabel([{'PFC Neurons'}; {'Organized by Peak F-Value'}]);
xlabel('Time Relative to Trial Start (s)');
hold on;
line([0 0], get(gca, 'ylim'), 'linestyle', '--', 'color', 'white', 'linewidth', 2);
% Odor Sorted by Previous Odor
erlyOdrPREVODRSortingMtx = [peakPrevOdrZ;1:length(peakOdrZ);normErlTrlOdrZ]';
erlyOdrPREVODRSortedMtx = sortrows(erlyOdrPREVODRSortingMtx);
erlyOdrPREVODRSortedIC = erlyOdrPREVODRSortedMtx(:,3:end);
subplot(3,3,6);
imagesc(pokeInIC(1).TimeBins,1:length(wholeTrialData),erlyOdrPREVODRSortedIC, [0 1]);
title('Current Odor Information Sorted by Previous Odor');
set(gca, 'xtick', 0, 'xticklabel', 'PokeIn', 'yticklabel', []);
ylabel([{'PFC Neurons'}; {'Organized by Peak F-Value'}]);
xlabel('Time Relative to Trial Start (s)');
hold on;
line([0 0], get(gca, 'ylim'), 'linestyle', '--', 'color', 'white', 'linewidth', 2);

% Previous Odor Sorted by Position
erlyPrevOdrPOSSortingMtx = [peakPosZ;1:length(peakPosZ);normErlTrlPrevOdrZ]';
erlyPrevOdrPOSSortedMtx = sortrows(erlyPrevOdrPOSSortingMtx);
erlyPrevOdrPOSSortedIC = erlyPrevOdrPOSSortedMtx(:,3:end);
subplot(3,3,7);
imagesc(pokeInIC(1).TimeBins,1:length(wholeTrialData),erlyPrevOdrPOSSortedIC, [0 1]);
title('Previous Odor Sorted by Position');
set(gca, 'xtick', 0, 'xticklabel', 'PokeIn', 'yticklabel', []);
ylabel([{'PFC Neurons'}; {'Organized by Peak F-Value'}]);
xlabel('Time Relative to Trial Start (s)');
hold on;
line([0 0], get(gca, 'ylim'), 'linestyle', '--', 'color', 'white', 'linewidth', 2);
% Previous Odor Sorted by CurrentOdor
erlyPrevOdrODRSortingMtx = [peakOdrZ;1:length(peakOdrZ);normErlTrlPrevOdrZ]';
erlyPrevOdrODRSortedMtx = sortrows(erlyPrevOdrODRSortingMtx);
erlyPrevOdrODRSortedIC = erlyPrevOdrODRSortedMtx(:,3:end);
subplot(3,3,8);
imagesc(pokeInIC(1).TimeBins,1:length(wholeTrialData),erlyPrevOdrODRSortedIC, [0 1]);
title('Previous Odor Sorted by Current Odor');
set(gca, 'xtick', 0, 'xticklabel', 'PokeIn', 'yticklabel', []);
ylabel([{'PFC Neurons'}; {'Organized by Peak F-Value'}]);
xlabel('Time Relative to Trial Start (s)');
hold on;
line([0 0], get(gca, 'ylim'), 'linestyle', '--', 'color', 'white', 'linewidth', 2);
% Previous Odor Sorted
erlyPrevOdrSortingMtx = [peakPrevOdrZ;1:length(peakOdrZ);normErlTrlPrevOdrZ]';
erlyPrevOdrSortedMtx = sortrows(erlyPrevOdrSortingMtx);
erlyPrevOdrSortedIC = erlyPrevOdrSortedMtx(:,3:end);
subplot(3,3,9);
imagesc(pokeInIC(1).TimeBins,1:length(wholeTrialData),erlyPrevOdrSortedIC, [0 1]);
title('Previous Odor Information');
set(gca, 'xtick', 0, 'xticklabel', 'PokeIn', 'yticklabel', []);
ylabel([{'PFC Neurons'}; {'Organized by Peak F-Value'}]);
xlabel('Time Relative to Trial Start (s)');
hold on;
line([0 0], get(gca, 'ylim'), 'linestyle', '--', 'color', 'white', 'linewidth', 2);


figure;
lateTrialPosZ = [pokeOutIC.PosZ];
lateTrialOdrZ = [pokeOutIC.OdorZ];
lateTrialPrevOdrZ = [pokeOutIC.PrevOdorZ];
normLtTrlPosZ = nan(size(lateTrialPosZ));
normLtTrlOdrZ = nan(size(lateTrialOdrZ));
normLtTrlPrevOdrZ = nan(size(lateTrialPrevOdrZ));
peakPosZ = nan(1,size(lateTrialPosZ,2));
peakOdrZ = nan(1,size(lateTrialOdrZ,2));
peakPrevOdrZ = nan(1,size(lateTrialPrevOdrZ,2));
for u = 1:size(lateTrialPosZ,2)
    curUniPosIC = lateTrialPosZ(:,u);
    maxPosIC = max(curUniPosIC);
    peakPosZ(u) = find(curUniPosIC==maxPosIC,1,'first');
    normLtTrlPosZ(:,u) = curUniPosIC./maxPosIC;
    
    curUniOdrIC = lateTrialOdrZ(:,u);
    maxOdrIC = max(curUniOdrIC);
    peakOdrZ(u) = find(curUniOdrIC==maxOdrIC,1,'first');
    normLtTrlOdrZ(:,u) = curUniOdrIC./maxOdrIC;
    
    curUniPrevOdrIC = lateTrialPrevOdrZ(:,u);
    maxPrevOdrIC = max(curUniPrevOdrIC);
    peakPrevOdrZ(u) = find(curUniPrevOdrIC==maxPrevOdrIC,1,'first');
    normLtTrlPrevOdrZ(:,u) = curUniPrevOdrIC./maxPrevOdrIC;
end

% Position Sorted by Position
latePosSortingMtx = [peakPosZ;1:length(peakPosZ);normLtTrlPosZ]';
latePosSortedMtx = sortrows(latePosSortingMtx);
latePosSortedIC = latePosSortedMtx(:,3:end);
subplot(3,3,1);
imagesc(pokeOutIC(1).TimeBins,1:length(wholeTrialData),latePosSortedIC, [0 1]);
title('Position Information');
set(gca, 'xtick', 0, 'xticklabel', 'PokeIn', 'yticklabel', []);
ylabel([{'PFC Neurons'}; {'Organized by Peak F-Value'}]);
xlabel('Time Relative to Trial End (s)');
hold on;
line([0 0], get(gca, 'ylim'), 'linestyle', '--', 'color', 'white', 'linewidth', 2);
% Position Sorted by Current Odor
latePosODRSortingMtx = [peakOdrZ;1:length(peakPosZ);normLtTrlPosZ]';
latePosODRSortedMtx = sortrows(latePosODRSortingMtx);
latePosODRSortedIC = latePosODRSortedMtx(:,3:end);
subplot(3,3,2);
imagesc(pokeOutIC(1).TimeBins,1:length(wholeTrialData),latePosODRSortedIC, [0 1]);
title('Position Sorted by Current Odor');
set(gca, 'xtick', 0, 'xticklabel', 'PokeIn', 'yticklabel', []);
ylabel([{'PFC Neurons'}; {'Organized by Peak F-Value'}]);
xlabel('Time Relative to Trial End (s)');
hold on;
line([0 0], get(gca, 'ylim'), 'linestyle', '--', 'color', 'white', 'linewidth', 2);
% Position Sorted by Previous Odor
latePosPREVODRSortingMtx = [peakPrevOdrZ;1:length(peakPrevOdrZ);normLtTrlPosZ]';
latePosPREVODRSortedMtx = sortrows(latePosPREVODRSortingMtx);
latePosPREVODRSortedIC = latePosPREVODRSortedMtx(:,3:end);
subplot(3,3,3);
imagesc(pokeOutIC(1).TimeBins,1:length(wholeTrialData),latePosPREVODRSortedIC, [0 1]);
title('Position Sorted by Previous Odor');
set(gca, 'xtick', 0, 'xticklabel', 'PokeIn', 'yticklabel', []);
ylabel([{'PFC Neurons'}; {'Organized by Peak F-Value'}]);
xlabel('Time Relative to Trial End (s)');
hold on;
line([0 0], get(gca, 'ylim'), 'linestyle', '--', 'color', 'white', 'linewidth', 2);

% Current Odor Sorted by Position
lateOdrPOSSortingMtx = [peakPosZ;1:length(peakPosZ);normLtTrlOdrZ]';
lateOdrPOSSortedMtx = sortrows(lateOdrPOSSortingMtx);
lateOdrPOSSortedIC = lateOdrPOSSortedMtx(:,3:end);
subplot(3,3,4);
imagesc(pokeOutIC(1).TimeBins,1:length(wholeTrialData),lateOdrPOSSortedIC, [0 1]);
title('Current Odor Sorted by Position');
set(gca, 'xtick', 0, 'xticklabel', 'PokeIn', 'yticklabel', []);
ylabel([{'PFC Neurons'}; {'Organized by Peak F-Value'}]);
xlabel('Time Relative to Trial End (s)');
hold on;
line([0 0], get(gca, 'ylim'), 'linestyle', '--', 'color', 'white', 'linewidth', 2);
% Odor Sorted by Odor
lateOdrSortingMtx = [peakOdrZ;1:length(peakOdrZ);normLtTrlOdrZ]';
lateOdrSortedMtx = sortrows(lateOdrSortingMtx);
lateOdrSortedIC = lateOdrSortedMtx(:,3:end);
subplot(3,3,5);
imagesc(pokeOutIC(1).TimeBins,1:length(wholeTrialData),lateOdrSortedIC, [0 1]);
title('Current Odor Information');
set(gca, 'xtick', 0, 'xticklabel', 'PokeIn', 'yticklabel', []);
ylabel([{'PFC Neurons'}; {'Organized by Peak F-Value'}]);
xlabel('Time Relative to Trial End (s)');
hold on;
line([0 0], get(gca, 'ylim'), 'linestyle', '--', 'color', 'white', 'linewidth', 2);
% Odor Sorted by Previous Odor
lateOdrPREVODRSortingMtx = [peakPrevOdrZ;1:length(peakOdrZ);normLtTrlOdrZ]';
lateOdrPREVODRSortedMtx = sortrows(lateOdrPREVODRSortingMtx);
lateOdrPREVODRSortedIC = lateOdrPREVODRSortedMtx(:,3:end);
subplot(3,3,6);
imagesc(pokeOutIC(1).TimeBins,1:length(wholeTrialData),lateOdrPREVODRSortedIC, [0 1]);
title('Current Odor Information Sorted by Previous Odor');
set(gca, 'xtick', 0, 'xticklabel', 'PokeIn', 'yticklabel', []);
ylabel([{'PFC Neurons'}; {'Organized by Peak F-Value'}]);
xlabel('Time Relative to Trial End (s)');
hold on;
line([0 0], get(gca, 'ylim'), 'linestyle', '--', 'color', 'white', 'linewidth', 2);

% Previous Odor Sorted by Position
latePrevOdrPOSSortingMtx = [peakPosZ;1:length(peakPosZ);normLtTrlPrevOdrZ]';
latePrevOdrPOSSortedMtx = sortrows(latePrevOdrPOSSortingMtx);
latePrevOdrPOSSortedIC = latePrevOdrPOSSortedMtx(:,3:end);
subplot(3,3,7);
imagesc(pokeOutIC(1).TimeBins,1:length(wholeTrialData),latePrevOdrPOSSortedIC, [0 1]);
title('Previous Odor Sorted by Position');
set(gca, 'xtick', 0, 'xticklabel', 'PokeIn', 'yticklabel', []);
ylabel([{'PFC Neurons'}; {'Organized by Peak F-Value'}]);
xlabel('Time Relative to Trial End (s)');
hold on;
line([0 0], get(gca, 'ylim'), 'linestyle', '--', 'color', 'white', 'linewidth', 2);
% Previous Odor Sorted by CurrentOdor
latePrevOdrODRSortingMtx = [peakOdrZ;1:length(peakOdrZ);normLtTrlPrevOdrZ]';
latePrevOdrODRSortedMtx = sortrows(latePrevOdrODRSortingMtx);
latePrevOdrODRSortedIC = latePrevOdrODRSortedMtx(:,3:end);
subplot(3,3,8);
imagesc(pokeOutIC(1).TimeBins,1:length(wholeTrialData),latePrevOdrODRSortedIC, [0 1]);
title('Previous Odor Sorted by Current Odor');
set(gca, 'xtick', 0, 'xticklabel', 'PokeIn', 'yticklabel', []);
ylabel([{'PFC Neurons'}; {'Organized by Peak F-Value'}]);
xlabel('Time Relative to Trial End (s)');
hold on;
line([0 0], get(gca, 'ylim'), 'linestyle', '--', 'color', 'white', 'linewidth', 2);
% Previous Odor Sorted
latePrevOdrSortingMtx = [peakPrevOdrZ;1:length(peakOdrZ);normLtTrlPrevOdrZ]';
latePrevOdrSortedMtx = sortrows(latePrevOdrSortingMtx);
latePrevOdrSortedIC = latePrevOdrSortedMtx(:,3:end);
subplot(3,3,9);
imagesc(pokeOutIC(1).TimeBins,1:length(wholeTrialData),latePrevOdrSortedIC, [0 1]);
title('Previous Odor Information');
set(gca, 'xtick', 0, 'xticklabel', 'PokeIn', 'yticklabel', []);
ylabel([{'PFC Neurons'}; {'Organized by Peak F-Value'}]);
xlabel('Time Relative to Trial End (s)');
hold on;
line([0 0], get(gca, 'ylim'), 'linestyle', '--', 'color', 'white', 'linewidth', 2);

%% Analysis #4 Position vs Odor
pokeInPosZ = [pokeInIC.PosZ];
pokeInOdrZ = [pokeInIC.OdorZ];

figure;
pokeInDiff = pokeInPosZ - pokeInOdrZ;

maxBiasNorm = nan(size(pokeInDiff));
maxBiasPos = nan(1,size(pokeInDiff,2));
maxBiasDir = nan(1,size(pokeInDiff,2));
maxOdrBiasNorm = nan(size(pokeInDiff));
maxOdrBiasPos = nan(1,size(pokeInDiff,2));
maxPosBiasNorm = nan(size(pokeInDiff));
maxPosBiasPos = nan(1,size(pokeInDiff,2));
for u = 1:length(ic)
    curUniBias = pokeInDiff(:,u);
    maxOdrBiasNorm(:,u) = curUniBias./abs(min(curUniBias));
    maxOdrBiasPos(u) = find(curUniBias==min(curUniBias),1,'first');
    maxPosBiasNorm(:,u) = curUniBias./max(curUniBias);
    maxPosBiasPos(u) = find(curUniBias==max(curUniBias),1,'first');
    if abs(min(curUniBias)) > max(curUniBias)
        maxBiasNorm(:,u) = curUniBias./abs(min(curUniBias));
        maxBiasPos(u) = find(curUniBias==min(curUniBias),1,'first');
        maxBiasDir(u) = -1;
    elseif abs(min(curUniBias)) < max(curUniBias)
        maxBiasNorm(:,u) = curUniBias./max(curUniBias);
        maxBiasPos(u) = find(curUniBias==max(curUniBias),1,'first');
        maxBiasDir(u) = 1;
    end
end
maxBiasDir(isnan(maxBiasDir)) = 0;
maxBiasPos(isnan(maxBiasDir)) = 1;
maxOdrBiasPos(isnan(maxOdrBiasPos)) = 1;
maxPosBiasPos(isnan(maxPosBiasPos)) = 1;

% erlyBiasSortingMtx = [peakFRsPI';1:length(maxBiasPos);maxBiasNorm]'; %% To Sort using peak pokeIn firing rate
erlyBiasSortingMtx = [maxBiasPos;1:length(maxBiasPos);maxBiasNorm]'; %% To sort using peak bias
erlyBiasSortedMtx = sortrows(erlyBiasSortingMtx);
erlyBiasSortedIC = erlyBiasSortedMtx(:,3:end);
a = subplot(5,1,1:4)
imagesc(pokeInIC(1).TimeBins,1:length(wholeTrialData),erlyBiasSortedIC, [-1 1]);
title('Information Bias');
set(gca, 'xtick', 0, 'xticklabel', 'PokeIn', 'yticklabel', []);
ylabel([{'PFC Neurons'}; {'Organized by Peak F-Value'}]);
xlabel('Time Relative to Trial Start (s)');
hold on;
line([0 0], get(gca, 'ylim'), 'linestyle', '--', 'color', 'white', 'linewidth', 2);

b = subplot(5,1,5)
plot(pokeInIC(1).TimeBins, mean(pokeInDiff'), 'color', 'k', 'linewidth', 2);
axis tight
hold on;
jbfill(pokeInIC(1).TimeBins', (mean(pokeInDiff')+ std(pokeInDiff')), (mean(pokeInDiff')-std(pokeInDiff')),[0 0 0], [0 0 0], 1, 0.5)
grid on;
box off
axis tight

figure;
erlyPosBiasSortingMtx = [maxPosBiasPos;1:length(maxOdrBiasPos);maxPosBiasNorm]';
erlyPosBiasSortedMtx = sortrows(erlyPosBiasSortingMtx);
erlyPosBiasSortedIC = erlyPosBiasSortedMtx(:,3:end);
subplot(2,2,1)
imagesc(pokeInIC(1).TimeBins,1:length(wholeTrialData),erlyPosBiasSortedIC, [-1 1]);
title('Position Information');
set(gca, 'xtick', -0.75:0.25:0.5, 'xticklabel', {'-0.75', '-0.5', '-0.25', 'PokeIn', '0.25', '0.5'}, 'yticklabel', []);
ylabel([{'PFC Neurons'}; {'Organized by Peak F-Value'}]);
xlabel('Time Relative to Trial Start (s)');
hold on;
line([0 0], get(gca, 'ylim'), 'linestyle', '--', 'color', 'white', 'linewidth', 2);

erlyOdrBiasSortingMtx = [maxOdrBiasPos;1:length(maxOdrBiasPos);maxOdrBiasNorm]';
erlyOdrBiasSortedMtx = sortrows(erlyOdrBiasSortingMtx);
erlyOdrBiasSortedIC = erlyOdrBiasSortedMtx(:,3:end);
subplot(2,2,2)
imagesc(pokeInIC(1).TimeBins,1:length(wholeTrialData),erlyOdrBiasSortedIC, [-1 1]);
title('Odor Information');
set(gca, 'xtick', -0.75:0.25:0.5, 'xticklabel', {'-0.75', '-0.5', '-0.25', 'PokeIn', '0.25', '0.5'}, 'yticklabel', []);
ylabel([{'PFC Neurons'}; {'Organized by Peak F-Value'}]);
xlabel('Time Relative to Trial Start (s)');
hold on;
line([0 0], get(gca, 'ylim'), 'linestyle', '--', 'color', 'white', 'linewidth', 2);

erlyPosBiasODRSortingMtx = [maxOdrBiasPos;1:length(maxOdrBiasPos);maxPosBiasNorm]';
erlyPosBiasODRSortedMtx = sortrows(erlyPosBiasODRSortingMtx);
erlyPosBiasODRSortedIC = erlyPosBiasODRSortedMtx(:,3:end);
subplot(2,2,3)
imagesc(pokeInIC(1).TimeBins,1:length(wholeTrialData),erlyPosBiasODRSortedIC, [-1 1]);
title('Position Information (sorted by Odor)');
set(gca, 'xtick', -0.75:0.25:0.5, 'xticklabel', {'-0.75', '-0.5', '-0.25', 'PokeIn', '0.25', '0.5'}, 'yticklabel', []);
ylabel([{'PFC Neurons'}; {'Organized by Peak F-Value'}]);
xlabel('Time Relative to Trial Start (s)');
hold on;
line([0 0], get(gca, 'ylim'), 'linestyle', '--', 'color', 'white', 'linewidth', 2);

erlyOdrBiasPOSSortingMtx = [maxPosBiasPos;1:length(maxOdrBiasPos);maxOdrBiasNorm]';
erlyOdrBiasPOSSortedMtx = sortrows(erlyOdrBiasPOSSortingMtx);
erlyOdrBiasPOSSortedIC = erlyOdrBiasPOSSortedMtx(:,3:end);
subplot(2,2,4)
imagesc(pokeInIC(1).TimeBins,1:length(wholeTrialData),erlyOdrBiasPOSSortedIC, [-1 1]);
title('Odor Information (sorted by Position)');
set(gca, 'xtick', -0.75:0.25:0.5, 'xticklabel', {'-0.75', '-0.5', '-0.25', 'PokeIn', '0.25', '0.5'}, 'yticklabel', []);
ylabel([{'PFC Neurons'}; {'Organized by Peak F-Value'}]);
xlabel('Time Relative to Trial Start (s)');
hold on;
line([0 0], get(gca, 'ylim'), 'linestyle', '--', 'color', 'white', 'linewidth', 2);

%
pokeOutPosZ = [pokeOutIC.PosZ];
pokeOutOdrZ = [pokeOutIC.OdorZ];

figure;
pokeOutDiff = pokeOutPosZ - pokeOutOdrZ;

maxBiasNorm = nan(size(pokeOutDiff));
maxBiasPos = nan(1,size(pokeOutDiff,2));
maxBiasDir = nan(1,size(pokeOutDiff,2));
maxOdrBiasNorm = nan(size(pokeOutDiff));
maxOdrBiasPos = nan(1,size(pokeOutDiff,2));
maxPosBiasNorm = nan(size(pokeOutDiff));
maxPosBiasPos = nan(1,size(pokeOutDiff,2));
for u = 1:length(ic)
    curUniBias = pokeOutDiff(:,u);
    maxOdrBiasNorm(:,u) = curUniBias./abs(min(curUniBias));
    maxOdrBiasPos(u) = find(curUniBias==min(curUniBias),1,'first');
    maxPosBiasNorm(:,u) = curUniBias./max(curUniBias);
    maxPosBiasPos(u) = find(curUniBias==max(curUniBias),1,'first');
    if abs(min(curUniBias)) > max(curUniBias)
        maxBiasNorm(:,u) = curUniBias./abs(min(curUniBias));
        maxBiasPos(u) = find(curUniBias==min(curUniBias),1,'first');
        maxBiasDir(u) = -1;
    elseif abs(min(curUniBias)) < max(curUniBias)
        maxBiasNorm(:,u) = curUniBias./max(curUniBias);
        maxBiasPos(u) = find(curUniBias==max(curUniBias),1,'first');
        maxBiasDir(u) = 1;
    end
end
maxBiasDir(isnan(maxBiasDir)) = 0;
maxBiasPos(isnan(maxBiasDir)) = 1;
maxOdrBiasPos(isnan(maxOdrBiasPos)) = 1;
maxPosBiasPos(isnan(maxPosBiasPos)) = 1;

erlyBiasSortingMtx = [maxBiasPos;1:length(maxBiasPos);maxBiasNorm]'; %% To sort using peak pokeOut bias latency
% erlyBiasSortingMtx = [peakFRsPO';1:length(maxBiasPos);maxBiasNorm]'; %% To sort using peak pokeOut firing rate
erlyBiasSortedMtx = sortrows(erlyBiasSortingMtx);
erlyBiasSortedIC = erlyBiasSortedMtx(:,3:end);
a = subplot(5,1,1:4)
imagesc(pokeOutIC(1).TimeBins,1:length(wholeTrialData),erlyBiasSortedIC, [-1 1]);
title('Information Bias');
set(gca, 'xtick', 0, 'xticklabel', 'pokeOut', 'yticklabel', []);
ylabel([{'PFC Neurons'}; {'Organized by Peak F-Value'}]);
xlabel('Time Relative to Trial Start (s)');
hold on;
line([0 0], get(gca, 'ylim'), 'linestyle', '--', 'color', 'white', 'linewidth', 2);

b = subplot(5,1,5)
plot(pokeOutIC(1).TimeBins, mean(pokeOutDiff'), 'color', 'k', 'linewidth', 2);
axis tight
hold on;
jbfill(pokeOutIC(1).TimeBins', (mean(pokeOutDiff')+ std(pokeOutDiff')), (mean(pokeOutDiff')-std(pokeOutDiff')),[0 0 0], [0 0 0], 1, 0.5)
grid on;
box off
axis tight

figure;
erlyPosBiasSortingMtx = [maxPosBiasPos;1:length(maxOdrBiasPos);maxPosBiasNorm]';
erlyPosBiasSortedMtx = sortrows(erlyPosBiasSortingMtx);
erlyPosBiasSortedIC = erlyPosBiasSortedMtx(:,3:end);
subplot(2,2,1)
imagesc(pokeOutIC(1).TimeBins,1:length(wholeTrialData),erlyPosBiasSortedIC, [-1 1]);
title('Position Information');
set(gca, 'xtick', -0.75:0.25:0.5, 'xticklabel', {'-0.75', '-0.5', '-0.25', 'pokeOut', '0.25', '0.5'}, 'yticklabel', []);
ylabel([{'PFC Neurons'}; {'Organized by Peak F-Value'}]);
xlabel('Time Relative to Trial Start (s)');
hold on;
line([0 0], get(gca, 'ylim'), 'linestyle', '--', 'color', 'white', 'linewidth', 2);

erlyOdrBiasSortingMtx = [maxOdrBiasPos;1:length(maxOdrBiasPos);maxOdrBiasNorm]';
erlyOdrBiasSortedMtx = sortrows(erlyOdrBiasSortingMtx);
erlyOdrBiasSortedIC = erlyOdrBiasSortedMtx(:,3:end);
subplot(2,2,2)
imagesc(pokeOutIC(1).TimeBins,1:length(wholeTrialData),erlyOdrBiasSortedIC, [-1 1]);
title('Odor Information');
set(gca, 'xtick', -0.75:0.25:0.5, 'xticklabel', {'-0.75', '-0.5', '-0.25', 'pokeOut', '0.25', '0.5'}, 'yticklabel', []);
ylabel([{'PFC Neurons'}; {'Organized by Peak F-Value'}]);
xlabel('Time Relative to Trial Start (s)');
hold on;
line([0 0], get(gca, 'ylim'), 'linestyle', '--', 'color', 'white', 'linewidth', 2);

erlyPosBiasODRSortingMtx = [maxOdrBiasPos;1:length(maxOdrBiasPos);maxPosBiasNorm]';
erlyPosBiasODRSortedMtx = sortrows(erlyPosBiasODRSortingMtx);
erlyPosBiasODRSortedIC = erlyPosBiasODRSortedMtx(:,3:end);
subplot(2,2,3)
imagesc(pokeOutIC(1).TimeBins,1:length(wholeTrialData),erlyPosBiasODRSortedIC, [-1 1]);
title('Position Information (sorted by Odor)');
set(gca, 'xtick', -0.75:0.25:0.5, 'xticklabel', {'-0.75', '-0.5', '-0.25', 'pokeOut', '0.25', '0.5'}, 'yticklabel', []);
ylabel([{'PFC Neurons'}; {'Organized by Peak F-Value'}]);
xlabel('Time Relative to Trial Start (s)');
hold on;
line([0 0], get(gca, 'ylim'), 'linestyle', '--', 'color', 'white', 'linewidth', 2);

erlyOdrBiasPOSSortingMtx = [maxPosBiasPos;1:length(maxOdrBiasPos);maxOdrBiasNorm]';
erlyOdrBiasPOSSortedMtx = sortrows(erlyOdrBiasPOSSortingMtx);
erlyOdrBiasPOSSortedIC = erlyOdrBiasPOSSortedMtx(:,3:end);
subplot(2,2,4)
imagesc(pokeOutIC(1).TimeBins,1:length(wholeTrialData),erlyOdrBiasPOSSortedIC, [-1 1]);
title('Odor Information (sorted by Position)');
set(gca, 'xtick', -0.75:0.25:0.5, 'xticklabel', {'-0.75', '-0.5', '-0.25', 'pokeOut', '0.25', '0.5'}, 'yticklabel', []);
ylabel([{'PFC Neurons'}; {'Organized by Peak F-Value'}]);
xlabel('Time Relative to Trial Start (s)');
hold on;
line([0 0], get(gca, 'ylim'), 'linestyle', '--', 'color', 'white', 'linewidth', 2);

%% Analysis #5 Peak FR vs Peak IC
pokeInPosZ = [pokeInIC.PosZ];
pokeInOdrZ = [pokeInIC.OdorZ];

posCorr = nan(1,size(pokeInPosZ,2));
odrCorr = nan(1,size(pokeInOdrZ,2));
for u = 1:size(pokeInPosZ,2)
    [r,~] = corrcoef(pokeInPosZ(:,u), normUniFRsPI(u,:));
    posCorr(u) = r(2);
    [r2,~] = corrcoef(pokeInOdrZ(:,u), normUniFRsPI(u,:));
    odrCorr(u) = r2(2);
end
figure;
subplot(2,1,1)
histogram(posCorr, -1:0.1:1)
subplot(2,1,2)
histogram(odrCorr, -1:0.1:1)
%% IC by Phase
% icPhase = [uniSumData.InformationContentSpikePhase];
% bands = fieldnames(icPhase);
% phaseBins = fieldnames(icPhase(1).(bands{1}));
% 
% phaseColors = [{'r'}, {'k'}, {'b'}];
% 
% for u = 1:length(icPhase)
%     f = figure;
%     spPlots = nan(1,length(bands));
%     for b = 1:length(bands)
%         spPlots(b) = subplot(length(bands),1,b);
%         hold on
%         switch b
%             case 1
%                 curBand = 'Theta';
%                 curOdrVals = icPhase(u).Theta;
%                 curPosVals = icPhase(u).Theta;
%             case 2
%                 curBand = 'Alpha';
%                 curOdrVals = icPhase(u).Alpha;
%                 curPosVals = icPhase(u).Alpha;
%             case 3
%                 curBand = 'Beta';
%                 curOdrVals = icPhase(u).Beta;
%                 curPosVals = icPhase(u).Beta;
%             case 4
%                 curBand = 'LowGamma';
%                 curOdrVals = icPhase(u).LowGamma;
%                 curPosVals = icPhase(u).LowGamma;
%             case 5
%                 curBand = 'HighGamma';
%                 curOdrVals = icPhase(u).HighGamma;
%                 curPosVals = icPhase(u).HighGamma;
%             case 6
%                 curBand = 'Ripple';
%                 curOdrVals = icPhase(u).Ripple;
%                 curPosVals = icPhase(u).Ripple;
%         end
%         for p = 1:length(phaseBins)
%             curOdrValBins = curOdrVals.(phaseBins{p}).EarlyTrial.PrevOdorZ;
%             curPosValBins = curPosVals.(phaseBins{p}).EarlyTrial.CurrPosZ;
% %             curOdrVals = icPhase(u).(bands{b}).(phaseBins{p}).EarlyTrial.PrevOdorZ;
% %             curPosVals = icPhase(p).(bands{b}).(phaseBins{p}).EarlyTrial.CurrPosZ;
%             plot(-0.8:0.001:0.5, curPosValBins - curOdrValBins, 'color', phaseColors{p});
% %             plot(-0.8:0.001:0.5, curPosVals, 'color', phaseColors{p});
% %             plot(-0.8:0.001:0.5, curOdrVals, 'color', phaseColors{p});
%         end
%         title([uniSumFiles{u} ' ' curBand], 'interpreter', 'none');
%         line(get(gca, 'xlim'), [0 0], 'color', 'k', 'linewidth', 1);
% %         set(gca, 'ylim', [-0.5 1]);
%         grid on;
%     end
%     linkaxes(spPlots, 'y');
%     
%     orient(gcf, 'tall');
%     orient(gcf, 'landscape');
% %     print
% %     print('-painters', gcf, '-dpdf', [uniSumFiles{u}(1:end-4) 'Phase']);
%     saveas(gcf, [uniSumFiles{u}(1:end-4) 'Phase.png'], 'png');
% end

%% IC Broken Down by Evoked Activity Periods & Modulation Periods
figure;
trialEpochStats = [uniSumData.TrialEpochStats];
ic = [uniSumData.InformationContent];
% Poke In Aligned
pokeIn = [ic.TrialPokeInSANSA];
piTimebins = pokeIn(1).TimeBins;
validPItimes = piTimebins<0.5;
piTimebins(~validPItimes) = [];
posPI = [pokeIn.PosZ];
posPI(~validPItimes,:) = [];
odrPI = [pokeIn.OdorZ];
odrPI(~validPItimes,:) = [];
prevOdrPI = [pokeIn.PrevOdorZ];
prevOdrPI(~validPItimes,:) = [];
pIn = subplot(1,2,1);
diffPlots = plot(piTimebins, posPI-odrPI, 'color', 'k', 'linewidth', 2);
title('PokeIn Aligned');
axis tight; grid on;
for a = 1:length(diffPlots)
    diffPlots(a).Color(4) = 0.25;
end
hold on;
plot(piTimebins,mean(posPI-prevOdrPI,2), 'color', 'r', 'linewidth', 3);
plot(piTimebins, mean(posPI-odrPI,2), 'color','g', 'linewidth', 2);
plot(piTimebins, mean(prevOdrPI-odrPI,2), 'color', 'c','linewidth', 2);
% Poke Out Aligned
pokeOut = [ic.TrialPokeOutSANSA];
poTimebins = pokeOut(1).TimeBins;
validPOtimes = poTimebins>-0.5;
poTimebins(~validPOtimes) = [];
posPO = [pokeOut.PosZ];
posPO(~validPOtimes,:) = [];
odrPO = [pokeOut.OdorZ];
odrPO(~validPOtimes,:) = [];
prevOdrPO = [pokeOut.PrevOdorZ];
prevOdrPO(~validPOtimes,:) = [];
pOut = subplot(1,2,2);
diffPlots = plot(poTimebins, posPO-odrPO, 'color', 'k', 'linewidth', 2);
title('PokeOut Aligned');
axis tight; grid on;
for a = 1:length(diffPlots)
    diffPlots(a).Color(4) = 0.25;
end
hold on;
plot(poTimebins,mean(posPO-prevOdrPO,2), 'color', 'r', 'linewidth', 3);
plot(poTimebins, mean(posPO-odrPO,2), 'color','g', 'linewidth', 2);
plot(poTimebins, mean(prevOdrPO-odrPO,2), 'color', 'c', 'linewidth', 2);

linkaxes([pIn, pOut], 'y');

%%
figure;
trialEpochStats = [uniSumData.TrialEpochStats];
ic = [uniSumData.InformationContent];
% Poke In Aligned
pokeIn = [ic.TrialPokeInSANSA];
piTimebins = pokeIn(1).TimeBins;
validPItimes = piTimebins<0.5;
piTimebins(~validPItimes) = [];
posPI = [pokeIn.PosZ];
posPI(~validPItimes,:) = [];
odrPI = [pokeIn.OdorZ];
odrPI(~validPItimes,:) = [];
prevOdrPI = [pokeIn.PrevOdorZ];
prevOdrPI(~validPItimes,:) = [];
pInPos = subplot(3,2,1);
posPlots = plot(piTimebins, posPI, 'color', 'k', 'linewidth', 2);
hold on;
plot(piTimebins, mean(posPI,2), 'color', 'r', 'linewidth', 2);
title('PokeIn Aligned');
ylabel('Position');
axis tight;
hold on;
pInOdr = subplot(3,2,3);
odrPlots = plot(piTimebins, odrPI, 'color', 'k', 'linewidth', 2);
hold on;
plot(piTimebins, mean(odrPI,2), 'color', 'r', 'linewidth', 2);
ylabel('Odor');
axis tight;
pInPrevOdr = subplot(3,2,5);
prevOdrPlots = plot(piTimebins, prevOdrPI, 'color', 'k', 'linewidth', 2);
hold on;
plot(piTimebins, mean(prevOdrPI,2), 'color', 'r', 'linewidth', 2);
ylabel('Previous Odor');
axis tight;
for a = 1:length(posPlots)
    posPlots(a).Color(4) = 0.25;
    odrPlots(a).Color(4) = 0.25;
    prevOdrPlots(a).Color(4) = 0.25;
end

% Poke Out Aligned
pokeOut = [ic.TrialPokeOutSANSA];
poTimebins = pokeOut(1).TimeBins;
validPOtimes = poTimebins>-0.5;
poTimebins(~validPOtimes) = [];
posPO = [pokeOut.PosZ];
posPO(~validPOtimes,:) = [];
odrPO = [pokeOut.OdorZ];
odrPO(~validPOtimes,:) = [];
prevOdrPO = [pokeOut.PrevOdorZ];
prevOdrPO(~validPOtimes,:) = [];
pOutPos = subplot(3,2,2);
posPlots = plot(poTimebins, posPO, 'color', 'k', 'linewidth', 2);
hold on;
plot(poTimebins, mean(posPO,2), 'color', 'r', 'linewidth', 2);
axis tight;
title('PokeOut Aligned');
axis tight;
hold on;
pOutOdr = subplot(3,2,4);
odrPlots = plot(poTimebins, odrPO, 'color', 'k', 'linewidth', 2);
hold on;
plot(poTimebins, mean(odrPO,2), 'color', 'r', 'linewidth', 2);
axis tight;
pOutPrevOdr = subplot(3,2,6);
prevOdrPlots = plot(poTimebins, prevOdrPO, 'color', 'k', 'linewidth', 2);
hold on;
plot(poTimebins, mean(prevOdrPO,2), 'color', 'r', 'linewidth', 2);
axis tight;
for a = 1:length(posPlots)
    posPlots(a).Color(4) = 0.25;
    odrPlots(a).Color(4) = 0.25;
    prevOdrPlots(a).Color(4) = 0.25;
end
linkaxes([pInPos,pInOdr,pInPrevOdr,pOutPos,pOutOdr,pOutPrevOdr], 'y');

figure;
piPlot = subplot(1,2,1);
piPosPlt = plot(piTimebins, mean(posPI,2), 'linewidth', 2, 'color', 'k');
hold on;
jbfill(piTimebins', (mean(posPI,2)+(std(posPI,1,2)./sqrt(size(posPI,2))))', (mean(posPI,2)-(std(posPI,1,2)./sqrt(size(posPI,2))))', [0 0 0], [0 0 0], 0, 0.5);
piOdrPlt = plot(piTimebins, mean(odrPI,2), 'linewidth', 2, 'color', 'r');
jbfill(piTimebins', (mean(odrPI,2)+(std(odrPI,1,2)./sqrt(size(odrPI,2))))', (mean(odrPI,2)-(std(odrPI,1,2)./sqrt(size(odrPI,2))))', [1 0 0], [1 0 0], 0, 0.5);
piPrevOdrPlt = plot(piTimebins, mean(prevOdrPI,2), 'linewidth', 2, 'color', 'c');
jbfill(piTimebins', (mean(prevOdrPI,2)+(std(prevOdrPI,1,2)./sqrt(size(prevOdrPI,2))))', (mean(prevOdrPI,2)-(std(prevOdrPI,1,2)./sqrt(size(prevOdrPI,2))))', [0 1 1], [0 1 1], 0, 0.5);
axis tight;
title('Information Content Relative to Poke In');
legend([piPosPlt, piOdrPlt, piPrevOdrPlt], {'Position', 'Odor', 'Previous Odor'}, 'location', 'northwest');
% legend([piPosPlt, piOdrPlt], {'Position', 'Odor'}, 'location', 'northwest');
ylabel({'Information Content'; 'Chance Normalized F-Ratios'});

poPlot = subplot(1,2,2);
plot(poTimebins, mean(posPO,2), 'linewidth', 2, 'color', 'k');
hold on;
jbfill(poTimebins', (mean(posPO,2)+(std(posPO,1,2)./sqrt(size(posPO,2))))', (mean(posPO,2)-(std(posPO,1,2)./sqrt(size(posPO,2))))', [0 0 0], [0 0 0], 0, 0.5);
plot(poTimebins, mean(odrPO,2), 'linewidth', 2, 'color', 'r');
jbfill(poTimebins', (mean(odrPO,2)+(std(odrPO,1,2)./sqrt(size(odrPO,2))))', (mean(odrPO,2)-(std(odrPO,1,2)./sqrt(size(odrPO,2))))', [1 0 0], [1 0 0], 0, 0.5);
plot(poTimebins, mean(prevOdrPO,2), 'linewidth', 2, 'color', 'c');
jbfill(poTimebins', (mean(prevOdrPO,2)+(std(prevOdrPO,1,2)./sqrt(size(prevOdrPO,2))))', (mean(prevOdrPO,2)-(std(prevOdrPO,1,2)./sqrt(size(prevOdrPO,2))))', [0 1 1], [0 1 1], 0, 0.5);
axis tight;
title('Information Content Relative to Poke Out');
linkaxes([piPlot, poPlot], 'y');
set(piPlot, 'ylim', [-1 5.5]);

%%
% 
% 
% trialPeriods = {'PreTrial', 'EarlyTrial', 'LateTrial', 'PostTrial','Reward','Error'};
% trlModMtxPos = nan(length(trialPeriods));
% trlModMtxOdr = nan(length(trialPeriods));
% for prd1 = 1:length(trialPeriods)
%     curPeriod = [trialEpochStats.(trialPeriods{prd1})];
%     curPrdOdor = cell2mat(cellfun(@(a){a'},{curPeriod.OdorSANSA}));
%     curPrdPos = cell2mat(cellfun(@(a){a'},{curPeriod.PositionSANSA}));
%     for prd2 = 1:length(trialPeriods)
%         curPrdLog = maxFRepoch==prd2;
%         trlModMtxOdr(prd2,prd1) = mean(curPrdOdor(2,curPrdLog)<newCritF)*75;
%         trlModMtxPos(prd2,prd1) = mean(curPrdPos(2,curPrdLog)<newCritF)*75;
%     end
% end