function IdentifySWRs_SM
%%
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
tetFileLog = cellfun(@(a)~isempty(a), regexp(fileNames, '_T([0-9]*)_SM.mat'));
tetFiles = fileNames(tetFileLog)';
behavFileLog = cellfun(@(a)~isempty(a), regexp(fileNames, '_BehaviorMatrix'));
load(fileNames{behavFileLog});

%% Create Smoothing Vector
w = gausswin(21);
w = w/sum(w);
%% Extract Raw Values & Compute RMS Power
ripBPF = cell(size(tetFiles))';
ripVolts = cell(size(tetFiles))';
ripRMS = cell(size(tetFiles))';
ripHilb = cell(size(tetFiles))';
ripTetIDs = cell(size(tetFiles));
for fl = 1:length(tetFiles)
    load(tetFiles{fl})
    samp = mode(diff(statMatrix(:,1)));
    wIndx = round((1/250)*(1/samp));
    fprintf('%s......', tetFiles{fl});
    ripCol = find(cellfun(@(a)~isempty(a), regexp(statMatrixColIDs, 'LFP_Ripple$')));
    ripVolts{fl} = statMatrix(:,2);
    ripRMS{fl} = conv(sqrt(conv(statMatrix(:,ripCol).^2, ones(wIndx,1)/wIndx, 'same')), w, 'same');
    ripBPF{fl} = statMatrix(:,ripCol);
    ripHilb{fl} = statMatrix(:,ripCol+1);
    ripTetIDs{fl} = statMatrixColIDs{ripCol};
    fprintf('done\n');
end

%% 
aggPower = mean(cell2mat(ripRMS),2);
bandPhase = cell2mat(ripHilb);
rmsThresh1 = (mean(aggPower) + (2*std(aggPower)));
rmsThresh2 = (mean(aggPower) + (4*std(aggPower)));
% aggPower = conv(aggPower, w, 'same');
abvThresh1 = aggPower>rmsThresh1;
abvThresh2 = aggPower>rmsThresh2;
epocWindows = [find(diff(abvThresh1)==1), find(diff(abvThresh1)==-1)];
dualThreshLog = false(size(epocWindows,1),1);
for epoc = 1:size(epocWindows,1)
    if sum(abvThresh2(epocWindows(epoc,1):epocWindows(epoc,2))) >=1
        dualThreshLog(epoc) = true;
    end
end
epocWindows(~dualThreshLog,:) = [];                                         % Comment out if not using dual thresholds
epocDur = diff(epocWindows,1,2);
interEpocInterval = epocWindows(2:end,1) - epocWindows(1:end-1,2);
% 
% Merge Ripples that Occur within 15ms of eachother
shortLatRipls = find(interEpocInterval<10);
for slr = 1:length(shortLatRipls)
    epocWindows(shortLatRipls(slr),2) = epocWindows(shortLatRipls(slr)+1,2);
end
epocWindows(shortLatRipls+1,:) = [];

epocDur = diff(epocWindows,1,2);
interEpocInterval = epocWindows(2:end,1) - epocWindows(1:end-1,2);

% epocWindows(epocDur<15,:) = [];                                             % Comment out if not using duration threshold

epocDur = diff(epocWindows,1,2);
interEpocInterval = epocWindows(2:end,1) - epocWindows(1:end-1,2);

figure;
subplot(2,2,1)
histogram(epocDur,0:1:100);
hold on;
line(repmat(mean(epocDur), [1,2]), get(gca, 'ylim'), 'color', 'r', 'linewidth', 2);
line(repmat(median(epocDur), [1,2]), get(gca, 'ylim'), 'color', 'r', 'linestyle', '--', 'linewidth', 2);
if mean(epocDur)>median(epocDur)
    set(gca, 'xtick', [median(epocDur), mean(epocDur)], 'xticklabel', [{'Median'}, {'Mean'}], 'xticklabelrotation', 90);
else
    set(gca, 'xtick', [mean(epocDur), median(epocDur)], 'xticklabel', [{'Mean'}, {'Median'}], 'xticklabelrotation', 90);
end
xlabel('Duration (ms)');
title([{'Ripple Duration'}; sprintf('Mean: %0.2f; Median: %0.2f', mean(epocDur), median(epocDur))]);
subplot(2,2,3)
histogram(interEpocInterval,0:10:500);
xlabel('Latency (ms)');
title('Inter-Ripple-Interval');
subplot(2,2,[2 4]);
corrScatPlot(epocDur(1:end-1), interEpocInterval, 'Duration', 'Latency',[]);

%% Format Output
tsVect = statMatrix(:,1);
swrIndices = tsVect(epocWindows);
swrVectLog = false(size(tsVect));
for rpl = 1:size(epocWindows,1)
    swrVectLog(epocWindows(rpl,1):epocWindows(rpl,2)) = true;
end
save('RippleEvents.mat', 'swrIndices', 'swrVectLog');
%%
chan = 2;

pokeCol = strcmp(behavMatrixColIDs, 'PokeEvents');
pokeEvents = [find(behavMatrix(:,pokeCol)==1), find(behavMatrix(:,pokeCol)==-1)];

scrsz = get(groot,'ScreenSize');
figure('position', [1 scrsz(4)*.45 scrsz(3) scrsz(4)*.47]); 
plot(ripRMS{chan}, 'linewidth', 2);
a = gca;
hold on;
plot(aggPower, 'linewidth', 2, 'color', 'g');
plot(a, ripBPF{chan}, 'color', 'r');
line(a,[0 length(ripRMS{chan})], repmat(mean(ripRMS{chan})+ std(ripRMS{chan})*2, [1 2]), 'linestyle', '-', 'color', 'black', 'linewidth', 2);
line(a,[0 length(ripRMS{chan})], repmat(mean(ripRMS{chan}) + std(ripRMS{chan})*4, [1 2]), 'linestyle', '--', 'color', 'black', 'linewidth', 2);
for pk = 1:length(pokeEvents)
    line(a, pokeEvents(pk,:), ones(1,2)*0.2, 'marker', '*', 'color', 'k');
end

for swr = 1:length(epocWindows)
    line(a,epocWindows(swr,:), ones(1,2)*0.25, 'marker', 'o', 'color', 'k');
end

figure('position', [1 1 scrsz(3) scrsz(4)/2]); 
plot(ripVolts{chan});
b = gca;
linkaxes([a b], 'x')

%% SOMETHING'S BROKE HERE
% peakLoc = nan(size(epocWindows,1),1);
% peakPwr = nan(size(epocWindows,1),1);
% phaseMtx = cell(size(epocWindows,1),1);
% phaseCorr = nan(size(epocWindows,1),1);
% parfor epoc = 1:size(epocWindows,1)
%     curEpocWindow = epocWindows(epoc,1):epocWindows(epoc,2);
%     [tmpPks, tmpLocs] = findpeaks(aggPower(curEpocWindow), 'MinPeakHeight', rmsThresh, 'MinPeakProminence', rmsThresh);
%     if isempty(tmpPks)
%         [tmpPks, tmpLocs] = findpeaks(aggPower(curEpocWindow));
%         %             if isempty(tmpPks)
%         %                 disp('EMPTY PEAKS')
%         %                 continue
%         % %                 tmpPks = aggPower(curEpocWindow);
%         % %                 tmpLocs = 1:length(tmpPks);
%         %             end
%     end
%     if ~isempty(tmpPks)
%         curEpocPhaseMtx = nan(size(bandPhase,2));
%         for tetA = 1:size(bandPhase,2)
%             for tetB = 1:size(bandPhase,2)
%                 curEpocPhaseMtx(tetA,tetB) = circ_corrcc(bandPhase(curEpocWindow,tetA), bandPhase(curEpocWindow,tetB));
%             end
%         end
%         phaseMtx{epoc} = curEpocPhaseMtx;
%         phaseCorr(epoc) = mean(curEpocPhaseMtx(triu(true(size(curEpocPhaseMtx)),1)));
%         peakPwr(epoc) = (max(tmpPks)-mean(aggPower))/std(aggPower); % z-norm peak power;
%         tmpEpocPkLoc = tmpLocs(tmpPks==max(tmpPks));
%         if size(tmpEpocPkLoc)>1
%             tmpEpocPkLoc = tmpEpocPkLoc(1);
%         end
%         tempEpocPkLoc =  epocWindows(epoc,1)+tmpEpocPkLoc;
%         peakLoc(epoc) = tempEpocPkLoc;
%     end
% end

% %% Edit the SWRs
% % Combine short SWRs
% 
% % Threshold SWRs with >=12.5ms duration 
% 
% %... may need to throw this into the loop above... check into this later
% %on.
% 
% %% 
% ripDur = diff(epocWindows,1,2)*samp;
% epocWindows(ripDur<0.0125,:) = [];
% peakPwr(ripDur<0.0125,:) = [];
% phaseCorr(ripDur<0.0125,:) = [];
% %%
% pokeCol = strcmp(behavMatrixColIDs, 'PokeEvents');
% 
% pokeEvents = [find(behavMatrix(:,pokeCol)==1), find(behavMatrix(:,pokeCol)==-1)];
% figure
% for pk = 1:length(pokeEvents)
%     line(pokeEvents(pk,:), ones(1,2), 'marker', '*', 'color', 'r');
%     hold on;
% end
% for swr = 1:length(epocWindows)
%     line(epocWindows(swr,:), ones(1,2)*2, 'marker', 'o', 'color', 'k');
% end
% 
% set(gca, 'ylim', [-1 4]);
% 
% %%
% figure;
% subplot(2,3,1)
% histogram(diff(epocWindows,1,2)*samp);
% 
% subplot(2,3,2)
% corrScatPlot(diff(epocWindows,1,2)*samp, peakPwr, 'SWR Dur', 'SWR Peak Power', []);
% 
% subplot(2,3,4)
% histogram(phaseCorr);
% 
% subplot(2,3,5)
% corrScatPlot(diff(epocWindows,1,2)*samp, phaseCorr, 'SWR Dur', 'Mean Coherence', []);
% 
% subplot(2,3,6)
% corrScatPlot(peakPwr, phaseCorr, 'SWR Peak Power', 'Mean Coherence', []);
% 
% 
% %%
% figure
% plot(repmat(statMatrix(:,1), [1, length(ripBPF)]), cell2mat(ripBPF))
% hold on;
% for epoc = 1:size(epocWindows,1)
%     line([statMatrix(epocWindows(epoc,1),1), statMatrix(epocWindows(epoc,2),1)], [0.2 0.2], 'marker', '*', 'color','k')
% end
% 


