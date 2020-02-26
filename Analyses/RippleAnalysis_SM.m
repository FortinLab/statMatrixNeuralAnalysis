% SWR Detection
%% Define Analysis Features
envProc = 'RMS';            % Enable for RMS determination of envelope
% envProc = 'HILB';           % Enable for abs(Hilbert) determinination of envelope
% %% Thresholds %% %
% Power Threshold %
powThresh1 = 0;             % Threshold 1
powThresh2 = 4;             % Threshold 2
% Duration Thresholds %
durThresh = 15;             % Duration Threshold
durThreshMrg = 15;          % Merge Duration Threshold
% Synchrony Threshold %
syncThresh = 0;             % Synchrony Threshold

% %% Parameters %% %
% Synchrony Duration %      % Not going to worry about this for now... as it will likely be computationally inefficient
syncWin = 10;
% Smoothing Vector %
w = gausswin(21);
w = w/sum(w);

% %% Trial Windows %% %
% Pre-Trial Window
preTrlWin = 500;
% Post-Trial Window
pstTrlWin = 700;
%% Load Data
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
tetFileLog = cellfun(@(a)~isempty(a), regexp(fileNames, '_T([0-9]*)_SM.mat'));
tetFiles = fileNames(tetFileLog)';
behavFileLog = cellfun(@(a)~isempty(a), regexp(fileNames, '_BehaviorMatrix'));
ensembleFileLog = cellfun(@(a)~isempty(a), regexp(fileNames, '_EnsembleMatrix'));
load(fileNames{behavFileLog});
load(fileNames{ensembleFileLog});
behavMatrixTrialStruct = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [0 0], 'PokeIn');

% poi = [behavMatrixTrialStruct.PokeOutIndex];
% ri = [behavMatrixTrialStruct.RewardIndex];
% if length(poi)~=length(ri)
%     ri(end:end+(length(poi)-length(ri)))=nan;
% end
% dif = ri-poi;
% perfLog = [behavMatrixTrialStruct.Performance]==1;
% isLog = [behavMatrixTrialStruct.TranspositionDistance]==0;
% iscLog = perfLog & isLog;
% isiLog = ~perfLog & isLog;
% oscLog = perfLog & ~isLog;
% osiLog = ~perfLog & ~isLog;
% 
% figure; 
% subplot(2,2,1); 
% histogram(dif(iscLog)); 
% title('InSeq Correct');
% subplot(2,2,2); 
% histogram(dif(isiLog)); 
% title('InSeq Incorrect');
% subplot(2,2,3); 
% histogram(dif(osiLog)); 
% title('OutSeq Incorrect');
% subplot(2,2,4); 
% histogram(dif(oscLog));
% title('OutSeq Correct');
% 
% annotation('textbox', 'position', [0.01 0.01 0.9 0.05], 'string',...
%     sprintf('%s', cd), 'linestyle', 'none', 'interpreter', 'none');
% 
% annotation('textbox', 'position', [0.1 0.95 0.8 0.05], 'string',...
%     '\bfReward Delivery Latency (relative to port withdrawal)', 'horizontalalignment', 'center',...
%     'linestyle', 'none', 'fontsize', 15);
%% Extract Raw Values & Compute RMS Power
ripBPF = nan(size(behavMatrix,1),size(tetFiles,1));
ripVolts = nan(size(behavMatrix,1),size(tetFiles,1));
ripRMS = nan(size(behavMatrix,1),size(tetFiles,1));
ripHilb = nan(size(behavMatrix,1),size(tetFiles,1));
ripTetIDs = cell(size(tetFiles));
for fl = 1:length(tetFiles)
    load(tetFiles{fl})
    samp = mode(diff(statMatrix(:,1)));
    wIndx = round((1/200)*(1/samp));
    fprintf('%s......', tetFiles{fl});
    ripCol = find(cellfun(@(a)~isempty(a), regexp(statMatrixColIDs, 'LFP_Ripple$')));
    ripVolts(:,fl) = statMatrix(:,2);
    if strcmp(envProc, 'RMS')
        ripRMS(:,fl) = conv(sqrt(conv(statMatrix(:,ripCol).^2, ones(wIndx,1)/wIndx, 'same')), w, 'same');
    elseif strcmp(envProc, 'HILB')
        ripRMS(:,fl) = conv(abs(hilbert(statMatrix(:,ripCol))),w,'same');
    end
    ripBPF(:,fl) = statMatrix(:,ripCol);
    ripHilb(:,fl) = statMatrix(:,ripCol+1);
    ripTetIDs{fl} = statMatrixColIDs{ripCol};
    fprintf('done\n');
end
%% Calculate Thresholds
% Aggregate Power
aggPower = mean(ripRMS,2);  aggType = 'Mean';               % Mean envelope
% aggPower = zscore(mean(ripRMS,2));  aggType = 'zMean';      % Z-Score Mean envelope
% aggPower = mean(zscore(ripRMS),2);  aggType = 'MeanZ';      % Mean Z-Score envelope

% Threshold Based on Mean +/- STD Aggregate Power
rmsThresh1 = (mean(aggPower) + (powThresh1*std(aggPower)));
rmsThresh2 = (mean(aggPower) + (powThresh2*std(aggPower)));

%% Identify Ripples
% Define putative ripple periods
abvThresh1 = aggPower>rmsThresh1;
epocWindows = [find(diff(abvThresh1)==1), find(diff(abvThresh1)==-1)];

% Apply secondary (peak power) threshold
abvThresh2 = aggPower>rmsThresh2;
dualThreshLog = false(size(epocWindows,1),1);
for epoc = 1:size(epocWindows,1)
    if sum(abvThresh2(epocWindows(epoc,1):epocWindows(epoc,2))) >=1
        dualThreshLog(epoc) = true;
    end
end
epocWindows(~dualThreshLog,:) = [];                                         % Comment out if not using dual thresholds

% Merge short latency ripples
interEpocInterval = epocWindows(2:end,1) - epocWindows(1:end-1,2);
slrLog = interEpocInterval<durThreshMrg;
shortLatRipls = find(slrLog);
mrgdNdxs = false(size(interEpocInterval));
for slr = 1:length(shortLatRipls)
    nxtNdx = shortLatRipls(slr)+find(slrLog(shortLatRipls(slr)+1:end)==0,1,'first');
    epocWindows(shortLatRipls(slr),2) = epocWindows(nxtNdx,2);
    mrgdNdxs(nxtNdx) = true;
end
epocWindows(mrgdNdxs,:) = [];

%% Quantify & Extract Ripples
% Determine the duration/timing of each event
epocDur = diff(epocWindows,1,2);

% Determine the Synchrony of each event
epocSync = nan(size(epocWindows,1),1);
% epocSyncConfMtx = cell(size(epocWindows,1),1);
epocRaw = cell(size(epocWindows,1),1);
epocBPF = cell(size(epocWindows,1),1);
epocSpike = cell(size(epocWindows,1),1);
for epoc = 1:size(epocWindows,1)
    epocRaw{epoc} = ripVolts(epocWindows(epoc,1):epocWindows(epoc,2),:);
    epocBPF{epoc} = ripBPF(epocWindows(epoc,1):epocWindows(epoc,2),:);
    epocSpike{epoc} = ensembleMatrix(epocWindows(epoc,1):epocWindows(epoc,2),2:end);
    tempConfMtx = nan(size(ripHilb,2));
    for e1 = 1:size(ripHilb,2)
        for e2 = 1:size(ripHilb,2)
            [tempConfMtx(e1,e2),~] = circ_corrcc(ripHilb(epocWindows(epoc,1):epocWindows(epoc,2),e1),...
                ripHilb(epocWindows(epoc,1):epocWindows(epoc,2),e2));
        end
    end
%     epocSyncConfMtx{epoc} = tempConfMtx;
    epocSync(epoc) = mean(tempConfMtx(triu(true(size(ripHilb,2)))));
end

%% Plot Duration vs Synchrony Correlation
figure;
dur = subplot(4,4,[1,5,9]);
histogram(epocDur, 30, 'Orientation', 'Horizontal');
set(gca, 'XDir', 'reverse');
line(get(gca, 'xlim'), repmat(mean(epocDur), [1,2]), 'color', 'r', 'linewidth', 2);
line(get(gca, 'xlim'), repmat(median(epocDur), [1,2]),'color', 'r', 'linestyle', '--', 'linewidth', 2);
if mean(epocDur)>median(epocDur)
    set(gca, 'ytick', [median(epocDur), mean(epocDur)], 'yticklabel', [{'Median'}, {'Mean'}]);
else
    set(gca, 'ytick', [mean(epocDur), median(epocDur)], 'yticklabel', [{'Mean'}, {'Median'}]);
end
xlabel('Count');
ylabel([{'Ripple Duration'}; sprintf('Mean: %0.2f; Median: %0.2f', mean(epocDur), median(epocDur))]);
corre = subplot(4,4,[2:4,6:8,10:12]);
corrScatPlot(epocSync, epocDur, [], [], []);
title('Session Wide Ripple Duration vs Mean Synchrony');
set(gca, 'xticklabel', [], 'yticklabel', []);
sync = subplot(4,4,14:16);
histogram(epocSync, 30);
xlabel('Synchrony (r)');
ylabel('Count');

linkaxes([dur, corre], 'y');
linkaxes([sync, corre], 'x');
annotation('textbox', 'position', [0.01 0.01 0.9 0.05], 'string',...
    sprintf('%s', cd), 'linestyle', 'none', 'interpreter', 'none');
annotation('textbox', 'position', [0.01 0.95 0.9 0.05], 'string',...
    sprintf('\\bfThreshold: \\rm+%i(+%i)SD; \\bfEnvelope = \\rm%s; \\bfAggregate = \\rm%s; \\bfMerge \\bfThreshold = \\rm%i ms', powThresh1 ,powThresh2, envProc, aggType, durThreshMrg),...
    'linestyle', 'none');

%% Extract Near-Trial Epocs
trialPokeTimes = [[behavMatrixTrialStruct.PokeInIndex]', [behavMatrixTrialStruct.PokeOutIndex]'];
trialRips = cell(size(trialPokeTimes,1),3);
trialRipLat = cell(size(trialPokeTimes,1),3);
trialRipsDur = cell(size(trialPokeTimes,1),3);
trialRipsSync = cell(size(trialPokeTimes,1),3);
for trl = 1:size(trialPokeTimes,1)
    preTrlLog = epocWindows(:,1)>(trialPokeTimes(trl,1)-preTrlWin) & epocWindows(:,1)<trialPokeTimes(trl,1);
    trialRips{trl,1} = epocWindows(preTrlLog,:);
    trialRipLat{trl,1} = epocWindows(preTrlLog,1) - trialPokeTimes(trl,1);
    trialRipsDur{trl,1} = epocDur(preTrlLog,:);
    trialRipsSync{trl,1} = epocSync(preTrlLog,:);
    
    trlLog = epocWindows(:,1)>trialPokeTimes(trl,1) & epocWindows(:,1)<trialPokeTimes(trl,2);
    trialRips{trl,2} = epocWindows(trlLog,:);
    trialRipLat{trl,2} = epocWindows(trlLog,1) - trialPokeTimes(trl,1);
    trialRipsDur{trl,2} = epocDur(trlLog,:);
    trialRipsSync{trl,2} = epocSync(trlLog,:);
    
    pstTrlLog = epocWindows(:,1)>trialPokeTimes(trl,2) & (epocWindows(:,1)<trialPokeTimes(trl,2)+pstTrlWin);
    trialRips{trl,3} = epocWindows(pstTrlLog,:);
    trialRipLat{trl,3} = epocWindows(pstTrlLog,1) - trialPokeTimes(trl,2);
    trialRipsDur{trl,3} = epocDur(pstTrlLog,:);
    trialRipsSync{trl,3} = epocSync(pstTrlLog,:);    
end

%% Plot Near-Trial Statistics
% Calculate Prcnt Trials w/Rip in Diff Period
figure;
prcntTrlsWrips = mean(cellfun(@(a)~isempty(a),trialRips));
subplot(3,1,1);
bar(prcntTrlsWrips);
set(gca, 'xticklabels', [], 'box', 'off');
ylabel('% Trials w/Ripple');
subplot(3,1,2);
BarPlotErrorbars([mean(cell2mat(trialRipsDur(:,1))),mean(cell2mat(trialRipsDur(:,2))),mean(cell2mat(trialRipsDur(:,3)))],...
    [std(cell2mat(trialRipsDur(:,1))),std(cell2mat(trialRipsDur(:,2))),std(cell2mat(trialRipsDur(:,3)))],[],[]);
set(gca, 'xticklabels', [], 'box', 'off');
ylabel('Duration (ms)');
subplot(3,1,3);
BarPlotErrorbars([mean(cell2mat(trialRipsSync(:,1))),mean(cell2mat(trialRipsSync(:,2))),mean(cell2mat(trialRipsSync(:,3)))],...
    [std(cell2mat(trialRipsSync(:,1))),std(cell2mat(trialRipsSync(:,2))),std(cell2mat(trialRipsSync(:,3)))],[],[]);
set(gca, 'xticklabels', [{'Pre-Trial'}, {'Trial'}, {'Post-Trial'}], 'box', 'off');
ylabel('Synchrony');
annotation('textbox', 'position', [0.01 0.01 0.9 0.05], 'string',...
    sprintf('%s', cd), 'linestyle', 'none', 'interpreter', 'none');
annotation('textbox', 'position', [0.01 0.95 0.9 0.05], 'string',...
    sprintf('\\bfThreshold: \\rm+%i(+%i)SD; \\bfEnvelope = \\rm%s; \\bfAggregate = \\rm%s; \\bfMerge \\bfThreshold = \\rm%i ms', powThresh1 ,powThresh2, envProc, aggType, durThreshMrg),...
    'linestyle', 'none');

%% Plot Break Down by InSeq/OutSeq & Perf
perfLog = [behavMatrixTrialStruct.Performance]==1;
isLog = [behavMatrixTrialStruct.TranspositionDistance]==0;
rips.Logicals.PerformanceLog = perfLog;
rips.Logicals.TransDist = [behavMatrixTrialStruct.TranspositionDistance];
rips.Logicals.OdorVect = [behavMatrixTrialStruct.Odor];
rips.Logicals.PositionVect = [behavMatrixTrialStruct.Position];
rips.TrialData.Events = trialRips;
rips.TrialData.Latency = trialRipLat;
rips.TrialData.Duration = trialRipsDur;
rips.TrialData.Synchrony = trialRipsSync;

figure;
iscLog = perfLog & isLog;
prcntTrlsWrips = mean(cellfun(@(a)~isempty(a),trialRips(iscLog,:)));
iscP = subplot(3,4,1);
bar(prcntTrlsWrips);
set(gca, 'xticklabels', [], 'box', 'off');
ylabel('% Trials w/Ripple');
title('InSeq Correct Trials');
iscD = subplot(3,4,5);
BarPlotErrorbars([mean(cell2mat(trialRipsDur(iscLog,1))),mean(cell2mat(trialRipsDur(iscLog,2))),mean(cell2mat(trialRipsDur(iscLog,3)))],...
    [std(cell2mat(trialRipsDur(iscLog,1))),std(cell2mat(trialRipsDur(iscLog,2))),std(cell2mat(trialRipsDur(iscLog,3)))],[],[]);
set(gca, 'xticklabels', [], 'box', 'off');
ylabel('Duration (ms)');
iscS = subplot(3,4,9);
BarPlotErrorbars([mean(cell2mat(trialRipsSync(iscLog,1))),mean(cell2mat(trialRipsSync(iscLog,2))),mean(cell2mat(trialRipsSync(iscLog,3)))],...
    [std(cell2mat(trialRipsSync(iscLog,1))),std(cell2mat(trialRipsSync(iscLog,2))),std(cell2mat(trialRipsSync(iscLog,3)))],[],[]);
set(gca, 'xticklabels', [{'Pre-Trial'}, {'Trial'}, {'Post-Trial'}], 'box', 'off');
ylabel('Synchrony');

isiLog = ~perfLog & isLog;
prcntTrlsWrips = mean(cellfun(@(a)~isempty(a),trialRips(isiLog,:)));
isiP = subplot(3,4,2);
bar(prcntTrlsWrips);
set(gca, 'xticklabels', [], 'box', 'off');
ylabel('% Trials w/Ripple');
title('InSeq Incorrect Trials');
isiD = subplot(3,4,6);
BarPlotErrorbars([mean(cell2mat(trialRipsDur(isiLog,1))),mean(cell2mat(trialRipsDur(isiLog,2))),mean(cell2mat(trialRipsDur(isiLog,3)))],...
    [std(cell2mat(trialRipsDur(isiLog,1))),std(cell2mat(trialRipsDur(isiLog,2))),std(cell2mat(trialRipsDur(isiLog,3)))],[],[]);
set(gca, 'xticklabels', [], 'box', 'off');
ylabel('Duration (ms)');
isiS = subplot(3,4,10);
BarPlotErrorbars([mean(cell2mat(trialRipsSync(isiLog,1))),mean(cell2mat(trialRipsSync(isiLog,2))),mean(cell2mat(trialRipsSync(isiLog,3)))],...
    [std(cell2mat(trialRipsSync(isiLog,1))),std(cell2mat(trialRipsSync(isiLog,2))),std(cell2mat(trialRipsSync(isiLog,3)))],[],[]);
set(gca, 'xticklabels', [{'Pre-Trial'}, {'Trial'}, {'Post-Trial'}], 'box', 'off');
ylabel('Synchrony');

oscLog = perfLog & ~isLog;
prcntTrlsWrips = mean(cellfun(@(a)~isempty(a),trialRips(oscLog,:)));
oscP = subplot(3,4,3);
bar(prcntTrlsWrips);
set(gca, 'xticklabels', [], 'box', 'off');
ylabel('% Trials w/Ripple');
title('OutSeq Correct Trials');
oscD = subplot(3,4,7);
BarPlotErrorbars([mean(cell2mat(trialRipsDur(oscLog,1))),mean(cell2mat(trialRipsDur(oscLog,2))),mean(cell2mat(trialRipsDur(oscLog,3)))],...
    [std(cell2mat(trialRipsDur(oscLog,1))),std(cell2mat(trialRipsDur(oscLog,2))),std(cell2mat(trialRipsDur(oscLog,3)))],[],[]);
set(gca, 'xticklabels', [], 'box', 'off');
ylabel('Duration (ms)');
oscS = subplot(3,4,11);
BarPlotErrorbars([mean(cell2mat(trialRipsSync(oscLog,1))),mean(cell2mat(trialRipsSync(oscLog,2))),mean(cell2mat(trialRipsSync(oscLog,3)))],...
    [std(cell2mat(trialRipsSync(oscLog,1))),std(cell2mat(trialRipsSync(oscLog,2))),std(cell2mat(trialRipsSync(oscLog,3)))],[],[]);
set(gca, 'xticklabels', [{'Pre-Trial'}, {'Trial'}, {'Post-Trial'}], 'box', 'off');
ylabel('Synchrony');

osiLog = ~perfLog & ~isLog;
prcntTrlsWrips = mean(cellfun(@(a)~isempty(a),trialRips(osiLog,:)));
osiP = subplot(3,4,4);
bar(prcntTrlsWrips);
set(gca, 'xticklabels', [], 'box', 'off');
ylabel('% Trials w/Ripple');
title('OutSeq Incorrect Trials');
osiD = subplot(3,4,8);
BarPlotErrorbars([mean(cell2mat(trialRipsDur(osiLog,1))),mean(cell2mat(trialRipsDur(osiLog,2))),mean(cell2mat(trialRipsDur(osiLog,3)))],...
    [std(cell2mat(trialRipsDur(osiLog,1))),std(cell2mat(trialRipsDur(osiLog,2))),std(cell2mat(trialRipsDur(osiLog,3)))],[],[]);
set(gca, 'xticklabels', [], 'box', 'off');
ylabel('Duration (ms)');
osiS = subplot(3,4,12);
BarPlotErrorbars([mean(cell2mat(trialRipsSync(osiLog,1))),mean(cell2mat(trialRipsSync(osiLog,2))),mean(cell2mat(trialRipsSync(osiLog,3)))],...
    [std(cell2mat(trialRipsSync(osiLog,1))),std(cell2mat(trialRipsSync(osiLog,2))),std(cell2mat(trialRipsSync(osiLog,3)))],[],[]);
set(gca, 'xticklabels', [{'Pre-Trial'}, {'Trial'}, {'Post-Trial'}], 'box', 'off');
ylabel('Synchrony');

annotation('textbox', 'position', [0.01 0.01 0.9 0.05], 'string',...
    sprintf('%s', cd), 'linestyle', 'none', 'interpreter', 'none');
annotation('textbox', 'position', [0.01 0.95 0.9 0.05], 'string',...
    sprintf('\\bfThreshold: \\rm+%i(+%i)SD; \\bfEnvelope = \\rm%s; \\bfAggregate = \\rm%s; \\bfMerge \\bfThreshold = \\rm%i ms', powThresh1 ,powThresh2, envProc, aggType, durThreshMrg),...
    'linestyle', 'none');

linkaxes([iscP isiP oscP osiP], 'xy');
linkaxes([iscD isiD oscD osiD], 'xy');
linkaxes([iscS isiS oscS osiS], 'xy');

%% Plot Breakdown of ISC by Odor
figure;
odorAlog = perfLog & isLog & [behavMatrixTrialStruct.Odor]==1;
prcntTrlsWrips = mean(cellfun(@(a)~isempty(a),trialRips(odorAlog,:)));
aP = subplot(3,5,1);
bar(prcntTrlsWrips);
set(gca, 'xticklabels', [], 'box', 'off');
ylabel('% Trials w/Ripple');
title('Odor A');
aD = subplot(3,5,6);
BarPlotErrorbars([mean(cell2mat(trialRipsDur(odorAlog,1))),mean(cell2mat(trialRipsDur(odorAlog,2))),mean(cell2mat(trialRipsDur(odorAlog,3)))],...
    [std(cell2mat(trialRipsDur(odorAlog,1))),std(cell2mat(trialRipsDur(odorAlog,2))),std(cell2mat(trialRipsDur(odorAlog,3)))],[],[]);
set(gca, 'xticklabels', [], 'box', 'off');
ylabel('Duration (ms)');
aS = subplot(3,5,11);
BarPlotErrorbars([mean(cell2mat(trialRipsSync(odorAlog,1))),mean(cell2mat(trialRipsSync(odorAlog,2))),mean(cell2mat(trialRipsSync(odorAlog,3)))],...
    [std(cell2mat(trialRipsSync(odorAlog,1))),std(cell2mat(trialRipsSync(odorAlog,2))),std(cell2mat(trialRipsSync(odorAlog,3)))],[],[]);
set(gca, 'xticklabels', [{'Pre-Trial'}, {'Trial'}, {'Post-Trial'}], 'box', 'off');
ylabel('Synchrony');
rips.A.RipWindows = trialRips(odorAlog,:);
rips.A.RipDur = trialRipsDur(odorAlog,:);
rips.A.RipSync = trialRipsSync(odorAlog,:);
rips.A.RipLat = trialRipLat(odorAlog,:);

odorBlog = perfLog & isLog & [behavMatrixTrialStruct.Odor]==2;
prcntTrlsWrips = mean(cellfun(@(a)~isempty(a),trialRips(odorBlog,:)));
bP = subplot(3,5,2);
bar(prcntTrlsWrips);
set(gca, 'xticklabels', [], 'box', 'off');
ylabel('% Trials w/Ripple');
title('Odor B');
bD = subplot(3,5,7);
BarPlotErrorbars([mean(cell2mat(trialRipsDur(odorBlog,1))),mean(cell2mat(trialRipsDur(odorBlog,2))),mean(cell2mat(trialRipsDur(odorBlog,3)))],...
    [std(cell2mat(trialRipsDur(odorBlog,1))),std(cell2mat(trialRipsDur(odorBlog,2))),std(cell2mat(trialRipsDur(odorBlog,3)))],[],[]);
set(gca, 'xticklabels', [], 'box', 'off');
ylabel('Duration (ms)');
bS = subplot(3,5,12);
BarPlotErrorbars([mean(cell2mat(trialRipsSync(odorBlog,1))),mean(cell2mat(trialRipsSync(odorBlog,2))),mean(cell2mat(trialRipsSync(odorBlog,3)))],...
    [std(cell2mat(trialRipsSync(odorBlog,1))),std(cell2mat(trialRipsSync(odorBlog,2))),std(cell2mat(trialRipsSync(odorBlog,3)))],[],[]);
set(gca, 'xticklabels', [{'Pre-Trial'}, {'Trial'}, {'Post-Trial'}], 'box', 'off');
ylabel('Synchrony');
rips.B.RipWindows = trialRips(odorBlog,:);
rips.B.RipDur = trialRipsDur(odorBlog,:);
rips.B.RipSync = trialRipsSync(odorBlog,:);
rips.B.RipLat = trialRipLat(odorBlog,:);

odorClog = perfLog & isLog & [behavMatrixTrialStruct.Odor]==3;
prcntTrlsWrips = mean(cellfun(@(a)~isempty(a),trialRips(odorClog,:)));
cP = subplot(3,5,3);
bar(prcntTrlsWrips);
set(gca, 'xticklabels', [], 'box', 'off');
ylabel('% Trials w/Ripple');
title('Odor C');
cD = subplot(3,5,8);
BarPlotErrorbars([mean(cell2mat(trialRipsDur(odorClog,1))),mean(cell2mat(trialRipsDur(odorClog,2))),mean(cell2mat(trialRipsDur(odorClog,3)))],...
    [std(cell2mat(trialRipsDur(odorClog,1))),std(cell2mat(trialRipsDur(odorClog,2))),std(cell2mat(trialRipsDur(odorClog,3)))],[],[]);
set(gca, 'xticklabels', [], 'box', 'off');
ylabel('Duration (ms)');
cS = subplot(3,5,13);
BarPlotErrorbars([mean(cell2mat(trialRipsSync(odorClog,1))),mean(cell2mat(trialRipsSync(odorClog,2))),mean(cell2mat(trialRipsSync(odorClog,3)))],...
    [std(cell2mat(trialRipsSync(odorClog,1))),std(cell2mat(trialRipsSync(odorClog,2))),std(cell2mat(trialRipsSync(odorClog,3)))],[],[]);
set(gca, 'xticklabels', [{'Pre-Trial'}, {'Trial'}, {'Post-Trial'}], 'box', 'off');
ylabel('Synchrony');
rips.C.RipWindows = trialRips(odorClog,:);
rips.C.RipDur = trialRipsDur(odorClog,:);
rips.C.RipSync = trialRipsSync(odorClog,:);
rips.C.RipLat = trialRipLat(odorClog,:);

odorDlog = perfLog & isLog & [behavMatrixTrialStruct.Odor]==4;
prcntTrlsWrips = mean(cellfun(@(a)~isempty(a),trialRips(odorDlog,:)));
dP = subplot(3,5,4);
bar(prcntTrlsWrips);
set(gca, 'xticklabels', [], 'box', 'off');
ylabel('% Trials w/Ripple');
title('Odor D');
dD = subplot(3,5,9);
BarPlotErrorbars([mean(cell2mat(trialRipsDur(odorDlog,1))),mean(cell2mat(trialRipsDur(odorDlog,2))),mean(cell2mat(trialRipsDur(odorDlog,3)))],...
    [std(cell2mat(trialRipsDur(odorDlog,1))),std(cell2mat(trialRipsDur(odorDlog,2))),std(cell2mat(trialRipsDur(odorDlog,3)))],[],[]);
set(gca, 'xticklabels', [], 'box', 'off');
ylabel('Duration (ms)');
dS = subplot(3,5,14);
BarPlotErrorbars([mean(cell2mat(trialRipsSync(odorDlog,1))),mean(cell2mat(trialRipsSync(odorDlog,2))),mean(cell2mat(trialRipsSync(odorDlog,3)))],...
    [std(cell2mat(trialRipsSync(odorDlog,1))),std(cell2mat(trialRipsSync(odorDlog,2))),std(cell2mat(trialRipsSync(odorDlog,3)))],[],[]);
set(gca, 'xticklabels', [{'Pre-Trial'}, {'Trial'}, {'Post-Trial'}], 'box', 'off');
ylabel('Synchrony');
rips.D.RipWindows = trialRips(odorDlog,:);
rips.D.RipDur = trialRipsDur(odorDlog,:);
rips.D.RipSync = trialRipsSync(odorDlog,:);
rips.D.RipLat = trialRipLat(odorDlog,:);

odorElog = perfLog & isLog & [behavMatrixTrialStruct.Odor]==5;
prcntTrlsWrips = mean(cellfun(@(a)~isempty(a),trialRips(odorElog,:)));
eP = subplot(3,5,5);
bar(prcntTrlsWrips);
set(gca, 'xticklabels', [], 'box', 'off');
ylabel('% Trials w/Ripple');
title('Odor E');
eD = subplot(3,5,10);
BarPlotErrorbars([mean(cell2mat(trialRipsDur(odorElog,1))),mean(cell2mat(trialRipsDur(odorElog,2))),mean(cell2mat(trialRipsDur(odorElog,3)))],...
    [std(cell2mat(trialRipsDur(odorElog,1))),std(cell2mat(trialRipsDur(odorElog,2))),std(cell2mat(trialRipsDur(odorElog,3)))],[],[]);
set(gca, 'xticklabels', [], 'box', 'off');
ylabel('Duration (ms)');
eS = subplot(3,5,15);
BarPlotErrorbars([mean(cell2mat(trialRipsSync(odorElog,1))),mean(cell2mat(trialRipsSync(odorElog,2))),mean(cell2mat(trialRipsSync(odorElog,3)))],...
    [std(cell2mat(trialRipsSync(odorElog,1))),std(cell2mat(trialRipsSync(odorElog,2))),std(cell2mat(trialRipsSync(odorElog,3)))],[],[]);
set(gca, 'xticklabels', [{'Pre-Trial'}, {'Trial'}, {'Post-Trial'}], 'box', 'off');
ylabel('Synchrony');
rips.E.RipWindows = trialRips(odorElog,:);
rips.E.RipDur = trialRipsDur(odorElog,:);
rips.E.RipSync = trialRipsSync(odorElog,:);
rips.E.RipLat = trialRipLat(odorElog,:);

annotation('textbox', 'position', [0.01 0.01 0.9 0.05], 'string',...
    sprintf('%s', cd), 'linestyle', 'none', 'interpreter', 'none');
annotation('textbox', 'position', [0.01 0.95 0.9 0.05], 'string',...
    sprintf('\\bfThreshold: \\rm+%i(+%i)SD; \\bfEnvelope = \\rm%s; \\bfAggregate = \\rm%s; \\bfMerge \\bfThreshold = \\rm%i ms', powThresh1 ,powThresh2, envProc, aggType, durThreshMrg),...
    'linestyle', 'none');

linkaxes([aP bP cP dP eP], 'xy');
linkaxes([aD bD cD dD eD], 'xy');
linkaxes([aS bS cS dS eS], 'xy');

%%
for e = 1:size(epocRaw,1)
    figure;
    
    spA = subplot(9,1,1:2);
    raw = plot(epocRaw{e}, 'color', 'k');
    for r = 1:length(raw)
        raw(r).Color(4) = 0.2;
    end
    title('Raw LFP')
    axis off
    spB = subplot(9,1,3:4);
    bpf = plot(epocBPF{e}, 'color', 'k');
    for b = 1:length(bpf)
        bpf(b).Color(4) = 0.2;
    end
    title('Bandpass LFP')
    axis off
    spC = subplot(9,1,5:8);
    [x, y] = find(epocSpike{e}~=0);
    scatter(x,y, '*k');
    title('Spiking')
    axis off
    
    spD = subplot(9,1,9);
    plot([ones(10,1); nan(size(epocRaw{e},1)-10,1)], 'k')
    text(3,1,'10ms', 'horizontalalignment', 'left', 'verticalalignment', 'bottom');
    axis off;
    
    linkaxes([spA, spB, spC, spD], 'x');
    
    annotation('textbox', 'position', [0.01 0.01 0.9 0.05], 'string',...
        sprintf('%s', cd), 'linestyle', 'none', 'interpreter', 'none');
    annotation('textbox', 'position', [0.1 0.01 0.8 0.05], 'string',...
        sprintf('%i / %i', e, size(epocRaw,1)), 'linestyle', 'none', 'interpreter', 'none', 'horizontalalignment', 'right');
    annotation('textbox', 'position', [0.01 0.95 0.9 0.05], 'string',...
        sprintf('\\bfThreshold: \\rm+%i(+%i)SD; \\bfEnvelope = \\rm%s; \\bfAggregate = \\rm%s; \\bfMerge \\bfThreshold = \\rm%i ms', powThresh1 ,powThresh2, envProc, aggType, durThreshMrg),...
        'linestyle', 'none');
    saveas(gcf, sprintf('Rips/Rip_%i(%i).bmp',e,size(epocRaw,1)));
    close gcf
end
