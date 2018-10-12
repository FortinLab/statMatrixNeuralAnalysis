function FvalAnal_PROTO
%%
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

%% Define Parameters
alignment = 'PokeIn';
windowStart = -1;
windowEnd = 1;
dataBinSize = 125;
%% Load Relevant Data
% 'trialInfo' : Matrix containing information about the identity and
%           outcome of each trial. Rows represent individual trials;
%           columns are organized thusly:
%               Column 1) Trial performance. 1 = Correct, 0 = Incorrect
%               Column 2) Trial InSeq log. 1 = InSeq, 0 = OutSeq
%               Column 3) Trial Position.
%               Column 4) Trial Odor.
[unitEpoch, unitIDs, lfpEpoch, lfpIDs, ~, eventTimeBins, trialInfo] = EpochExtraction_SM(alignment, windowStart, windowEnd, 'org', 'TiUTr', 'lfpBand', 'Beta');

%% Bin the unit data
samp = mode(diff(eventTimeBins));
binLims = linspace(windowStart, windowEnd, (windowEnd-windowStart)/(dataBinSize*samp)+2);
binnedUnitEpoch = nan(length(binLims)-1, size(unitEpoch,2), size(unitEpoch,3));
for uni = 1:length(unitIDs)
    for trl = 1:size(trialInfo,1)
        binnedUnitEpoch(:,uni,trl) = histcounts(eventTimeBins(unitEpoch(:,uni,trl)==1), binLims);
    end
end
binnedTimeBins = binLims(2:end)-(dataBinSize*samp/2);

%% Analysis #1 (Previous trial Odor Vs. Upcoming Trial Position)
corrTrlPos = nan(size(trialInfo,1),1);
inCorrTrlPos = nan(size(trialInfo,1),1);
corrTrlOdr = nan(size(trialInfo,1),1);
inCorrTrlOdr = nan(size(trialInfo,1),1);
isFtrOStrlPos = nan(size(trialInfo,1),1);
isFtrOStrlOdr = nan(size(trialInfo,1),1);
for trl = 2:size(trialInfo,1)
    if trialInfo(trl,1)==1 && (trialInfo(trl,3)>1 && trialInfo(trl,3)<length(unique(trialInfo(:,3)))) && (trialInfo(trl,3) - trialInfo(trl-1,3) == 1) % && trialInfo(trl-1,2)~=1
        corrTrlPos(trl) = trialInfo(trl,3);
        corrTrlOdr(trl) = trialInfo(trl-1,4);
    elseif trialInfo(trl,1)~=1 && (trialInfo(trl,3)>1 && trialInfo(trl,3)<length(unique(trialInfo(:,3)))) && (trialInfo(trl,3) - trialInfo(trl-1,3) == 1)
        inCorrTrlPos(trl) = trialInfo(trl,3);
        inCorrTrlOdr(trl) = trialInfo(trl-1,4);
    end
    if trialInfo(trl-1,1)==1 && trialInfo(trl-1,2)==0 && trialInfo(trl,1)==1 && (trialInfo(trl,3)-trialInfo(trl-1,3)==1)
        isFtrOStrlPos(trl) = trialInfo(trl,3);
        isFtrOStrlOdr(trl) = trialInfo(trl-1,4);
    end
end
corrTrlEnsmbl = binnedUnitEpoch(:,:,~isnan(corrTrlPos));
corrTrlPos = corrTrlPos(~isnan(corrTrlPos));
corrTrlOdr = corrTrlOdr(~isnan(corrTrlOdr));            

ftrOStrlEnsmbl = binnedUnitEpoch(:,:,~isnan(isFtrOStrlPos));
isFtrOStrlPos = isFtrOStrlPos(~isnan(isFtrOStrlPos));
isFtrOStrlOdr = isFtrOStrlOdr(~isnan(isFtrOStrlOdr));

unitPosFvals = nan(length(binnedTimeBins), length(unitIDs));
unitOdrFvals = nan(length(binnedTimeBins), length(unitIDs));
ftrOStrlPosFvals = nan(length(binnedTimeBins), length(unitIDs));
ftrOStrlOdrFvals = nan(length(binnedTimeBins), length(unitIDs));
for u = 1:length(unitIDs)
    for t = 1:length(binnedTimeBins)        
        [~,tablePos,~] = anova1(reshape(corrTrlEnsmbl(t,u,:), [size(corrTrlEnsmbl,3),1]), corrTrlPos, 'off');
        if ~isempty(tablePos{2,5})
            unitPosFvals(t,u) = tablePos{2,5};
        else
            unitPosFvals(t,u) = 1;
        end
        
        [~,ftrOSpos,~] = anova1(reshape(ftrOStrlEnsmbl(t,u,:), [size(ftrOStrlEnsmbl,3),1]), isFtrOStrlPos, 'off');
        if ~isempty(ftrOSpos{2,5})
            ftrOStrlPosFvals(t,u) = ftrOSpos{2,5};
        else
            ftrOStrlPosFvals(t,u) = 1;
        end
        
        [~,tableOdr,~] = anova1(reshape(corrTrlEnsmbl(t,u,:), [size(corrTrlEnsmbl,3),1]), corrTrlOdr, 'off');
        if ~isempty(tableOdr{2,5})
            unitOdrFvals(t,u) = tableOdr{2,5};
        else
            unitOdrFvals(t,u) = 1;
        end
        
        [~,ftrOSodr,~] = anova1(reshape(ftrOStrlEnsmbl(t,u,:), [size(ftrOStrlEnsmbl,3),1]), isFtrOStrlOdr, 'off');
        if ~isempty(ftrOSodr{2,5})
            ftrOStrlOdrFvals(t,u) = ftrOSodr{2,5};
        else
            ftrOStrlOdrFvals(t,u) = 1;
        end
    end
end
unitPosFvals(isnan(unitPosFvals)) = 1;
ftrOStrlPosFvals(isnan(ftrOStrlPosFvals)) = 1;
unitOdrFvals(isnan(unitOdrFvals)) = 1;
ftrOStrlOdrFvals(isnan(ftrOStrlOdrFvals)) = 1;

% % % % % % % % Maybe remove eventually
unitPosFvals(unitPosFvals>10) = 10;
ftrOStrlPosFvals(ftrOStrlPosFvals>10) = 10;
unitOdrFvals(unitOdrFvals>10) = 10;
ftrOStrlOdrFvals(ftrOStrlOdrFvals>10) = 10;

figure
subplot(2,1,1)
plot(binnedTimeBins, median(unitPosFvals - unitOdrFvals,2), 'color', 'k');
hold on
plot(binnedTimeBins, ((std(unitPosFvals - unitOdrFvals,1,2)./(length(unitIDs)-1))+median(unitPosFvals - unitOdrFvals,2)), 'linestyle', ':', 'color', 'black');
plot(binnedTimeBins, median(unitPosFvals - unitOdrFvals,2)-(std(unitPosFvals - unitOdrFvals,1,2)./(length(unitIDs)-1)), 'linestyle', ':', 'color', 'black');
line(get(gca, 'xlim'), [0 0], 'color', 'k', 'linewidth', 1);
set(gca, 'xtick', binLims, 'ylim', [-2 2]);
grid on;
title('Information Bias: Current Position vs Previous Odor (All Trials)', 'interpreter', 'none')
ylabel({'Fval Difference'; 'Positive Values = Position Bias'; 'Negative Values = Previous Odor Bias'});
xlabel(sprintf('Time relative to %s(s)', alignment));

subplot(2,1,2)
plot(binnedTimeBins, median(ftrOStrlPosFvals - ftrOStrlOdrFvals,2), 'color', 'k');
hold on;
plot(binnedTimeBins, ((std(ftrOStrlPosFvals - ftrOStrlOdrFvals,1,2)./(length(unitIDs)-1))+median(ftrOStrlPosFvals - ftrOStrlOdrFvals,2)), 'linestyle', ':', 'color', 'black');
plot(binnedTimeBins, median(ftrOStrlPosFvals - ftrOStrlOdrFvals,2)-(std(ftrOStrlPosFvals - ftrOStrlOdrFvals,1,2)./(length(unitIDs)-1)), 'linestyle', ':', 'color', 'black');
line(get(gca, 'xlim'), [0 0], 'color', 'k', 'linewidth', 1);
set(gca, 'xtick', binLims, 'ylim', [-2 2]);
grid on;
title('Information Bias: Trial After OutSeq Only');
ylabel({'Fval Difference'; 'Positive Values = Position Bias'; 'Negative Values = Previous Odor Bias'});
xlabel(sprintf('Time relative to %s(s)', alignment));

% figure
% subplot(3,1,1);
% plot(binnedTimeBins, unitPosFvals);
% subplot(3,1,2);
% plot(binnedTimeBins, unitOdrFvals);
% subplot(3,1,3);
% plot(binnedTimeBins, unitPosFvals-unitOdrFvals);
% legend(unitIDs)
%% Analysis #2 (Current Odor Vs. Current Position)
perfLog = trialInfo(:,1)==1;
aLog = (trialInfo(:,3)==1 | trialInfo(:,4)==1);

% Select data
corrTrlEnsmbl= binnedUnitEpoch(:,:,perfLog);
corrTrlPos = trialInfo(perfLog,3);
corrTrlOdr = trialInfo(perfLog,4);
% 
corrTrlEnsmblSansA = binnedUnitEpoch(:,:,perfLog & ~aLog);
corrTrlPosSansA = trialInfo(perfLog & ~aLog,3);
corrTrlOdrSansA = trialInfo(perfLog & ~aLog,4);

unitPosFvals = nan(length(binnedTimeBins), length(unitIDs));
unitOdrFvals = nan(length(binnedTimeBins), length(unitIDs));
unitPosFvalsSansA = nan(length(binnedTimeBins), length(unitIDs));
unitOdrFvalsSansA = nan(length(binnedTimeBins), length(unitIDs));
for u = 1:length(unitIDs)
    for t = 1:length(binnedTimeBins)        
        [~,tablePos,~] = anova1(reshape(corrTrlEnsmbl(t,u,:), [size(corrTrlEnsmbl,3),1]), corrTrlPos, 'off');
        unitPosFvals(t,u) = tablePos{2,5};
        [~,tablePosSansA,~] = anova1(reshape(corrTrlEnsmblSansA(t,u,:), [size(corrTrlEnsmblSansA,3),1]), corrTrlPosSansA, 'off');
        unitPosFvalsSansA(t,u) = tablePosSansA{2,5};
        
        [~,tableOdr,~] = anova1(reshape(corrTrlEnsmbl(t,u,:), [size(corrTrlEnsmbl,3),1]), corrTrlOdr, 'off');
        unitOdrFvals(t,u) = tableOdr{2,5};
        [~,tableOdrSansA,~] = anova1(reshape(corrTrlEnsmblSansA(t,u,:), [size(corrTrlEnsmblSansA,3),1]), corrTrlOdrSansA, 'off');
        unitOdrFvalsSansA(t,u) = tableOdrSansA{2,5};
    end
end
unitPosFvals(isnan(unitPosFvals)) = 1;
unitOdrFvals(isnan(unitOdrFvals)) = 1;
unitPosFvalsSansA(isnan(unitPosFvalsSansA)) = 1;
unitOdrFvalsSansA(isnan(unitOdrFvalsSansA)) = 1;

% % % % % % % % Maybe remove eventually
unitPosFvals(unitPosFvals>10) = 10;
unitPosFvalsSansA(unitPosFvalsSansA>10) = 10;
unitOdrFvals(unitOdrFvals>10) = 10;
unitOdrFvalsSansA(unitOdrFvalsSansA>10) = 10;

figure
subplot(2,1,1)
plot(binnedTimeBins, mean(unitPosFvals - unitOdrFvals,2), 'color','k');
hold on
plot(binnedTimeBins, ((std(unitPosFvals - unitOdrFvals,1,2)./(length(unitIDs)-1))+mean(unitPosFvals - unitOdrFvals,2)), 'linestyle', ':', 'color', 'black');
plot(binnedTimeBins, mean(unitPosFvals - unitOdrFvals,2)-(std(unitPosFvals - unitOdrFvals,1,2)./(length(unitIDs)-1)), 'linestyle', ':', 'color', 'black');
line(get(gca, 'xlim'), [0 0], 'color', 'k', 'linewidth', 1);
set(gca, 'xtick', binLims, 'ylim', [-2 2]);
grid on;
title('Information Bias: Current Position vs Current Odor (All Trials; Mean +/- SEM)', 'interpreter', 'none')
ylabel({'Fval Difference'; 'Positive Values = Position Bias'; 'Negative Values = Previous Odor Bias'});
xlabel(sprintf('Time relative to %s(s)', alignment));
subplot(2,1,2)
plot(binnedTimeBins, mean(unitPosFvalsSansA - unitOdrFvalsSansA,2), 'color','k');
hold on
plot(binnedTimeBins, ((std(unitPosFvalsSansA - unitOdrFvalsSansA,1,2)./(length(unitIDs)-1))+mean(unitPosFvalsSansA - unitOdrFvalsSansA,2)), 'linestyle', ':', 'color', 'black');
plot(binnedTimeBins, mean(unitPosFvalsSansA - unitOdrFvalsSansA,2)-(std(unitPosFvalsSansA - unitOdrFvalsSansA,1,2)./(length(unitIDs)-1)), 'linestyle', ':', 'color', 'black');
line(get(gca, 'xlim'), [0 0], 'color', 'k', 'linewidth', 1);
set(gca, 'xtick', binLims, 'ylim', [-2 2]);
grid on;
title('Information Bias: Current Position vs Current Odor (Odors B-D; Mean +/- SEM)', 'interpreter', 'none')
ylabel({'Fval Difference'; 'Positive Values = Position Bias'; 'Negative Values = Previous Odor Bias'});
xlabel(sprintf('Time relative to %s(s)', alignment));
% 
% figure
% subplot(3,1,1);
% plot(binnedTimeBins, unitPosFvals);
% subplot(3,1,2);
% plot(binnedTimeBins, unitOdrFvals);
% subplot(3,1,3);
% plot(binnedTimeBins, unitPosFvals-unitOdrFvals);
% legend(unitIDs)

%% Analysis #3 (Current Odor vs. Next Position)
corrTrlPos = nan(size(trialInfo,1),1);
corrTrlOdr = nan(size(trialInfo,1),1);
corrTrlPosSansA = nan(size(trialInfo,1),1);
corrTrlOdrSansA = nan(size(trialInfo,1),1);
for trl = 1:size(trialInfo,1)-1
    if trialInfo(trl,1)==1 && trialInfo(trl+1,1)==1 && (trialInfo(trl+1,3) - trialInfo(trl,3) == 1)
        corrTrlPos(trl) = trialInfo(trl+1,3);
        corrTrlOdr(trl) = trialInfo(trl,4);
    end
    if trialInfo(trl,1)==1 && trialInfo(trl+1,1)==1 && (trialInfo(trl+1,3) - trialInfo(trl,3) == 1) && trialInfo(trl,3)~=1
        corrTrlPosSansA(trl) = trialInfo(trl+1,3);
        corrTrlOdrSansA(trl) = trialInfo(trl,4);
    end
end

% Select data
corrTrlEnsmbl = binnedUnitEpoch(:,:,~isnan(corrTrlPos));
corrTrlPos = trialInfo(~isnan(corrTrlPos),3);
corrTrlOdr = trialInfo(~isnan(corrTrlOdr),4);
% 
corrTrlEnsmblSansA = binnedUnitEpoch(:,:,~isnan(corrTrlPosSansA));
corrTrlPosSansA = trialInfo(~isnan(corrTrlPosSansA),3);
corrTrlOdrSansA = trialInfo(~isnan(corrTrlOdrSansA),4);

unitPosFvals = nan(length(binnedTimeBins), length(unitIDs));
unitOdrFvals = nan(length(binnedTimeBins), length(unitIDs));
unitPosFvalsSansA = nan(length(binnedTimeBins), length(unitIDs));
unitOdrFvalsSansA = nan(length(binnedTimeBins), length(unitIDs));
for u = 1:length(unitIDs)
    for t = 1:length(binnedTimeBins)        
        [~,tablePos,~] = anova1(reshape(corrTrlEnsmbl(t,u,:), [size(corrTrlEnsmbl,3),1]), corrTrlPos, 'off');
        unitPosFvals(t,u) = tablePos{2,5};
        [~,tablePosSansA,~] = anova1(reshape(corrTrlEnsmblSansA(t,u,:), [size(corrTrlEnsmblSansA,3),1]), corrTrlPosSansA, 'off');
        unitPosFvalsSansA(t,u) = tablePosSansA{2,5};
        
        [~,tableOdr,~] = anova1(reshape(corrTrlEnsmbl(t,u,:), [size(corrTrlEnsmbl,3),1]), corrTrlOdr, 'off');
        unitOdrFvals(t,u) = tableOdr{2,5};
        [~,tableOdrSansA,~] = anova1(reshape(corrTrlEnsmblSansA(t,u,:), [size(corrTrlEnsmblSansA,3),1]), corrTrlOdrSansA, 'off');
        unitOdrFvalsSansA(t,u) = tableOdrSansA{2,5};
    end
end
unitPosFvals(isnan(unitPosFvals)) = 1;
unitOdrFvals(isnan(unitOdrFvals)) = 1;
unitPosFvalsSansA(isnan(unitPosFvalsSansA)) = 1;
unitOdrFvalsSansA(isnan(unitOdrFvalsSansA)) = 1;

% % % % % % % % Maybe remove eventually
unitPosFvals(unitPosFvals>10) = 10;
unitPosFvalsSansA(unitPosFvalsSansA>10) = 10;
unitOdrFvals(unitOdrFvals>10) = 10;
unitOdrFvalsSansA(unitOdrFvalsSansA>10) = 10;

figure
subplot(2,1,1)
plot(binnedTimeBins, mean(unitPosFvals - unitOdrFvals,2), 'color', 'k');
hold on
plot(binnedTimeBins, ((std(unitPosFvals - unitOdrFvals,1,2)./(length(unitIDs)-1))+mean(unitPosFvals - unitOdrFvals,2)), 'linestyle', ':', 'color', 'black');
plot(binnedTimeBins, mean(unitPosFvals - unitOdrFvals,2)-(std(unitPosFvals - unitOdrFvals,1,2)./(length(unitIDs)-1)), 'linestyle', ':', 'color', 'black');
line(get(gca, 'xlim'), [0 0], 'color', 'k', 'linewidth', 1);
set(gca, 'xtick', binLims, 'ylim', [-2 2]);
grid on;
title('Information Bias: Next Position vs Current Odor (All Trials; Mean +/- SEM)', 'interpreter', 'none')
ylabel({'Fval Difference'; 'Positive Values = Position Bias'; 'Negative Values = Previous Odor Bias'});
xlabel(sprintf('Time relative to %s(s)', alignment));
subplot(2,1,2)
plot(binnedTimeBins, mean(unitPosFvalsSansA - unitOdrFvalsSansA,2), 'color', 'k');
hold on
plot(binnedTimeBins, ((std(unitPosFvalsSansA - unitOdrFvalsSansA,1,2)./(length(unitIDs)-1))+mean(unitPosFvalsSansA - unitOdrFvalsSansA,2)), 'linestyle', ':', 'color', 'black');
plot(binnedTimeBins, mean(unitPosFvalsSansA - unitOdrFvalsSansA,2)-(std(unitPosFvalsSansA - unitOdrFvalsSansA,1,2)./(length(unitIDs)-1)), 'linestyle', ':', 'color', 'black');
line(get(gca, 'xlim'), [0 0], 'color', 'k', 'linewidth', 1);
set(gca, 'xtick', binLims, 'ylim', [-2 2]);
grid on;
title('Information Bias: Next Position vs Current Odor (Odors B&C; Mean +/- SEM)', 'interpreter', 'none')
ylabel({'Fval Difference'; 'Positive Values = Position Bias'; 'Negative Values = Previous Odor Bias'});
xlabel(sprintf('Time relative to %s(s)', alignment));

%% Evaluate Beta Power Modulation
% Compute Trialwise RMS
betaPower = nan(length(binnedTimeBins), length(lfpIDs), size(trialInfo,1));
for tet = 1:length(lfpIDs)
    for trl = 1:size(trialInfo,1)
        curRMS = conv(lfpEpoch(:,tet,trl).^2, ones(floor((1/30)/samp),1), 'same');
        for bin = 1:length(binLims)-1
            timeLog = eventTimeBins>=binLims(bin) & eventTimeBins<binLims(bin+1);
            betaPower(bin,tet,trl) = mean(curRMS(timeLog));
        end
    end
end

lfpTetIDs = cellfun(@(b)b{1},cellfun(@(a)strsplit(a,'_'), lfpIDs, 'uniformoutput', 0), 'uniformoutput', 0);
uniTetIDs = cellfun(@(b)b{1},cellfun(@(a)strsplit(a,'-'), unitIDs, 'uniformoutput',0), 'uniformoutput', 0);

for uni = 1:length(unitIDs)
    curUni = uniTetIDs(uni);
    lfpLog = strcmp(curUni, lfpTetIDs);
    curUniBinnedSpiking = binnedUnitEpoch(:,uni,:);
    curTetBinnedBetaPower = betaPower(:,lfpLog,:);
    
    figure;
    subplot(2,1,1)
    corrScatPlot(curUniBinnedSpiking(:)./(dataBinSize*samp), curTetBinnedBetaPower(:), 'Firing Rate', 'Beta Power', []);
    title(unitIDs{uni});
    subplot(2,1,2) 
    plot(binnedTimeBins,unitPosFvals(:,uni));
    hold on
    plot(binnedTimeBins,unitOdrFvals(:,uni));
    plot(binnedTimeBins,unitPosFvals(:,uni) - unitOdrFvals(:,uni));
    title(unitIDs{uni});
    legend('Pos', 'Odr', 'Diff')
end

%% OLD CODE
%%
%%
%% Compute Trialwise RMS
tetRMSvals = nan(size(trialInfo,1), length(lfpIDs));
for tet = 1:length(lfpIDs)
    for trl = 1:size(trialInfo,1)
        tetRMSvals(trl,tet) = rms(lfpEpoch(:,tet,trl));
    end
end
avgTrlRMS = median(tetRMSvals,2);
       
highBetaLog = avgTrlRMS>median(avgTrlRMS);
