function FvalAnalSlide_PROTO
origDir = cd;
[fileDir] = uigetdir(origDir);
if fileDir==0
    disp('Analysis Cancelled')
    return
else
    cd(fileDir)
end

%% Load Relevant Data
% 'trialInfo' : Matrix containing information about the identity and
%           outcome of each trial. Rows represent individual trials;
%           columns are organized thusly:
%               Column 1) Trial performance. 1 = Correct, 0 = Incorrect
%               Column 2) Trial InSeq log. 1 = InSeq, 0 = OutSeq
%               Column 3) Trial Position.
%               Column 4) Trial Odor.
%               Column 5) Poke Duration
alignment = 'PokeIn';
windowStart = -0.9;
windowEnd = 0.6;
[unitEpoch, unitIDs, lfpEpoch, lfpIDs, ~, eventTimeBins, trialInfo] = EpochExtraction_SM(alignment, windowStart, windowEnd, 'org', 'TiUTr', 'lfpBand', 'All', 'lfpData', 'Phase');

%% Define Parameters
dataBinSize = 200;
timePeriod = eventTimeBins(1+dataBinSize/2:length(eventTimeBins)-(dataBinSize/2));

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
corrTrlEnsmbl = unitEpoch(:,:,~isnan(corrTrlPos));
corrTrlPos = corrTrlPos(~isnan(corrTrlPos));
corrTrlOdr = corrTrlOdr(~isnan(corrTrlOdr));            

ftrOStrlEnsmbl = unitEpoch(:,:,~isnan(isFtrOStrlPos));
isFtrOStrlPos = isFtrOStrlPos(~isnan(isFtrOStrlPos));
isFtrOStrlOdr = isFtrOStrlOdr(~isnan(isFtrOStrlOdr));

unitPosFvals = nan(length(eventTimeBins)-dataBinSize, length(unitIDs));
unitOdrFvals = nan(length(eventTimeBins)-dataBinSize, length(unitIDs));
ftrOStrlPosFvals = nan(length(eventTimeBins)-dataBinSize, length(unitIDs));
ftrOStrlOdrFvals = nan(length(eventTimeBins)-dataBinSize, length(unitIDs));
for u = 1:length(unitIDs)
    tic
    fprintf('Running %s...', unitIDs{u});
    parfor t = 1+dataBinSize/2:length(eventTimeBins)-(dataBinSize/2)
        timeMask = false(length(eventTimeBins),1);
        timeMask(t-(dataBinSize/2):t+dataBinSize/2) = true;
        [~,tablePos,~] = anova1(reshape(sum(corrTrlEnsmbl(timeMask,u,:)), [size(corrTrlEnsmbl,3),1]), corrTrlPos, 'off');
        if ~isempty(tablePos{2,5})
            unitPosFvals(t,u) = tablePos{2,5};
        else
            unitPosFvals(t,u) = 1;
        end
        
        [~,ftrOSpos,~] = anova1(reshape(sum(ftrOStrlEnsmbl(timeMask,u,:)), [size(ftrOStrlEnsmbl,3),1]), isFtrOStrlPos, 'off');
        if ~isempty(ftrOSpos{2,5})
            ftrOStrlPosFvals(t,u) = ftrOSpos{2,5};
        else
            ftrOStrlPosFvals(t,u) = 1;
        end
        
        [~,tableOdr,~] = anova1(reshape(sum(corrTrlEnsmbl(timeMask,u,:)), [size(corrTrlEnsmbl,3),1]), corrTrlOdr, 'off');
        if ~isempty(tableOdr{2,5})
            unitOdrFvals(t,u) = tableOdr{2,5};
        else
            unitOdrFvals(t,u) = 1;
        end
        
        [~,ftrOSodr,~] = anova1(reshape(sum(ftrOStrlEnsmbl(timeMask,u,:)), [size(ftrOStrlEnsmbl,3),1]), isFtrOStrlOdr, 'off');
        if ~isempty(ftrOSodr{2,5})
            ftrOStrlOdrFvals(t,u) = ftrOSodr{2,5};
        else
            ftrOStrlOdrFvals(t,u) = 1;
        end
    end
%     fprintf('Completed in ');
    toc
end
% Clean up the data (remove blank entries at start & replace nans with 1s)
unitPosFvals(1:dataBinSize/2,:) = [];
unitPosFvals(isnan(unitPosFvals)) = 1;

ftrOStrlPosFvals(1:dataBinSize/2,:) = [];
ftrOStrlPosFvals(isnan(ftrOStrlPosFvals)) = 1;

unitOdrFvals(1:dataBinSize/2,:) = [];
unitOdrFvals(isnan(unitOdrFvals)) = 1;

ftrOStrlOdrFvals(1:dataBinSize/2,:) = [];
ftrOStrlOdrFvals(isnan(ftrOStrlOdrFvals)) = 1;

% Plot Everything
figure
subplot(2,1,1)
plot(timePeriod, mean(unitPosFvals - unitOdrFvals,2), 'color', 'k');
hold on
plot(timePeriod, ((std(unitPosFvals - unitOdrFvals,1,2)./(length(unitIDs)-1))+mean(unitPosFvals - unitOdrFvals,2)), 'linestyle', ':', 'color', 'black');
plot(timePeriod, mean(unitPosFvals - unitOdrFvals,2)-(std(unitPosFvals - unitOdrFvals,1,2)./(length(unitIDs)-1)), 'linestyle', ':', 'color', 'black');
line(get(gca, 'xlim'), [0 0], 'color', 'k', 'linewidth', 1);
set(gca, 'ylim', [-2 2]);
grid on;
title(sprintf('Information Bias: Current Position vs Previous Odor (All Trials) (%ims window)', dataBinSize), 'interpreter', 'none')
ylabel({'Fval Difference'; 'Positive Values = Position Bias'; 'Negative Values = Previous Odor Bias'});
xlabel(sprintf('Time relative to %s(s)', alignment));

subplot(2,1,2)
plot(timePeriod, mean(ftrOStrlPosFvals - ftrOStrlOdrFvals,2), 'color', 'k');
hold on;
plot(timePeriod, ((std(ftrOStrlPosFvals - ftrOStrlOdrFvals,1,2)./(length(unitIDs)-1))+mean(ftrOStrlPosFvals - ftrOStrlOdrFvals,2)), 'linestyle', ':', 'color', 'black');
plot(timePeriod, mean(ftrOStrlPosFvals - ftrOStrlOdrFvals,2)-(std(ftrOStrlPosFvals - ftrOStrlOdrFvals,1,2)./(length(unitIDs)-1)), 'linestyle', ':', 'color', 'black');
line(get(gca, 'xlim'), [0 0], 'color', 'k', 'linewidth', 1);
set(gca, 'ylim', [-2 2]);
grid on;
title('Information Bias: Trial After OutSeq Only');
ylabel({'Fval Difference'; 'Positive Values = Position Bias'; 'Negative Values = Previous Odor Bias'});
xlabel(sprintf('Time relative to %s(s)', alignment));
orient(gcf, 'tall');
orient(gcf, 'landscape');
annotation(gcf,'textbox', [0.01 0.01 0.96 0.03], 'FitBoxToText','off', 'string', cd, 'interpreter', 'none', 'linestyle', 'none');
drawnow
 
% figure
% subplot(3,1,1);
% plot(timePeriod, unitPosFvals);
% subplot(3,1,2);
% plot(timePeriod, unitOdrFvals);
% subplot(3,1,3);
% plot(timePeriod, unitPosFvals-unitOdrFvals);
% legend(unitIDs)


%% Analysis #2 Previous trial Odor Vs. Upcoming Trial Position: Skips vs Repeats
transInSkpPOS = nan(size(trialInfo,1),1);
transInSkpODR = nan(size(trialInfo,1),1);
transInRepPOS = nan(size(trialInfo,1),1);
transInRepODR = nan(size(trialInfo,1),1);
for trl = 2:size(trialInfo,1)
    if trialInfo(trl,1)==1 && trialInfo(trl-1,1)==1 && (trialInfo(trl,3)-trialInfo(trl-1,3))==1 && trialInfo(trl-1,2)~=1
        if (trialInfo(trl-1,3)-trialInfo(trl-1,4)) >= 1
            transInRepPOS(trl) = trialInfo(trl,3);
            transInRepODR(trl) = trialInfo(trl-1,4);
        elseif (trialInfo(trl-1,3)-trialInfo(trl-1,4)) <= -1
            transInSkpPOS(trl) = trialInfo(trl,3);
            transInSkpODR(trl) = trialInfo(trl-1,4);
        end
    elseif trialInfo(trl,3)~=1 && trialInfo(trl,1)==1 && trialInfo(trl-1,1)==1 && (trialInfo(trl,3)-trialInfo(trl-3,3))==1 && trialInfo(trl-1,2)==1
        transInRepPOS(trl) = trialInfo(trl,3);
        transInRepODR(trl) = trialInfo(trl-1,4);
        transInSkpPOS(trl) = trialInfo(trl,3);
        transInSkpODR(trl) = trialInfo(trl-1,4);
    end
end
skipTrlEnsmbl = unitEpoch(:,:,~isnan(transInSkpPOS));
skipTrlPos = transInSkpPOS(~isnan(transInSkpPOS));
skipTrlOdr = transInSkpODR(~isnan(transInSkpODR));

repTrlEnsmbl = unitEpoch(:,:,~isnan(transInRepPOS));
repTrlPos = transInRepPOS(~isnan(transInRepPOS));
repTrlOdr = transInRepODR(~isnan(transInRepODR));

skipPosFvals = nan(length(eventTimeBins)-dataBinSize, length(unitIDs));
skipOdrFvals = nan(length(eventTimeBins)-dataBinSize, length(unitIDs));
repPosFvals = nan(length(eventTimeBins)-dataBinSize, length(unitIDs));
repOdrFvals = nan(length(eventTimeBins)-dataBinSize, length(unitIDs));
for u = 1:length(unitIDs)
    tic
    fprintf('Running %s...', unitIDs{u});
    for t = 1+dataBinSize/2:length(eventTimeBins)-(dataBinSize/2)
        timeMask = false(length(eventTimeBins),1);
        timeMask(t-(dataBinSize/2):t+dataBinSize/2) = true;
        [~,skipPos,~] = anova1(reshape(sum(skipTrlEnsmbl(timeMask,u,:)), [size(skipTrlEnsmbl,3),1]), skipTrlPos, 'off');
        if ~isempty(skipPos{2,5})
            skipPosFvals(t,u) = skipPos{2,5};
        else
            skipPosFvals(t,u) = 1;
        end
        
        [~,skipOdr,~] = anova1(reshape(sum(skipTrlEnsmbl(timeMask,u,:)), [size(skipTrlEnsmbl,3),1]), skipTrlOdr, 'off');
        if ~isempty(skipOdr{2,5})
            skipOdrFvals(t,u) = skipOdr{2,5};
        else
            skipOdrFvals(t,u) = 1;
        end
        
        [~,repPos,~] = anova1(reshape(sum(repTrlEnsmbl(timeMask,u,:)), [size(repTrlEnsmbl,3),1]), repTrlPos, 'off');
        if ~isempty(repPos{2,5})
            repPosFvals(t,u) = repPos{2,5};
        else
            repPosFvals(t,u) = 1;
        end
        
        [~,repOdr,~] = anova1(reshape(sum(repTrlEnsmbl(timeMask,u,:)), [size(repTrlEnsmbl,3),1]), repTrlOdr, 'off');
        if ~isempty(repOdr{2,5})
            repOdrFvals(t,u) = repOdr{2,5};
        else
            repOdrFvals(t,u) = 1;
        end
    end
%     fprintf('Completed in ');
    toc
end
% Clean up the data (remove blank entries at start & replace nans with 1s)
skipPosFvals(1:dataBinSize/2,:) = [];
skipPosFvals(isnan(skipPosFvals)) = 1;

skipOdrFvals(1:dataBinSize/2,:) = [];
skipOdrFvals(isnan(skipOdrFvals)) = 1;

repPosFvals(1:dataBinSize/2,:) = [];
repPosFvals(isnan(repPosFvals)) = 1;

repOdrFvals(1:dataBinSize/2,:) = [];
repOdrFvals(isnan(repOdrFvals)) = 1;

% Plot Everything
figure
subplot(2,1,1)
plot(timePeriod, mean(skipPosFvals - skipOdrFvals,2), 'color', 'k');
hold on
plot(timePeriod, ((std(skipPosFvals - skipOdrFvals,1,2)./(length(unitIDs)-1))+mean(skipPosFvals - skipOdrFvals,2)), 'linestyle', ':', 'color', 'black');
plot(timePeriod, mean(skipPosFvals - skipOdrFvals,2)-(std(skipPosFvals - skipOdrFvals,1,2)./(length(unitIDs)-1)), 'linestyle', ':', 'color', 'black');
line(get(gca, 'xlim'), [0 0], 'color', 'k', 'linewidth', 1);
set(gca, 'ylim', [-2 2]);
grid on;
title(sprintf('Information Bias: Current Position vs Previous Odor (Skips Only) (%ims window)', dataBinSize), 'interpreter', 'none')
ylabel({'Fval Difference'; 'Positive Values = Position Bias'; 'Negative Values = Previous Odor Bias'});
xlabel(sprintf('Time relative to %s(s)', alignment));

subplot(2,1,2)
plot(timePeriod, mean(repPosFvals - repOdrFvals,2), 'color', 'k');
hold on;
plot(timePeriod, ((std(repPosFvals - repOdrFvals,1,2)./(length(unitIDs)-1))+mean(repPosFvals - repOdrFvals,2)), 'linestyle', ':', 'color', 'black');
plot(timePeriod, mean(repPosFvals - repOdrFvals,2)-(std(repPosFvals - repOdrFvals,1,2)./(length(unitIDs)-1)), 'linestyle', ':', 'color', 'black');
line(get(gca, 'xlim'), [0 0], 'color', 'k', 'linewidth', 1);
set(gca, 'ylim', [-2 2]);
grid on;
title('Information Bias: Current Position vs Previous Odor (Repeats Only)');
ylabel({'Fval Difference'; 'Positive Values = Position Bias'; 'Negative Values = Previous Odor Bias'});
xlabel(sprintf('Time relative to %s(s)', alignment));
orient(gcf, 'tall');
orient(gcf, 'landscape');
annotation(gcf,'textbox', [0.01 0.01 0.96 0.03], 'FitBoxToText','off', 'string', cd, 'interpreter', 'none', 'linestyle', 'none');
drawnow
 
% figure
% subplot(3,1,1);
% plot(timePeriod, unitPosFvals);
% subplot(3,1,2);
% plot(timePeriod, unitOdrFvals);
% subplot(3,1,3);
% plot(timePeriod, unitPosFvals-unitOdrFvals);
% legend(unitIDs)
%% Analysis #3: Previous Odor vs Upcoming Position: Organized by Phase
% Parse lfpIDs
lfpIDparts = cellfun(@(b)[b(1);b(3)], cellfun(@(a)strsplit(a, '_'), lfpIDs, 'uniformoutput', 0), 'uniformoutput', 0);
lfpIDparts = [lfpIDparts{:}];
bands = unique(lfpIDparts(2,:));
bands(strcmp(bands, 'Raw')) = [];
phaseBins = linspace(-pi,pi,8);

phaseColors = [{'r'}, {'k'}, {'b'}];
phaseBinLabel = cell(1,length(phaseBins)-1);
for p = 1:length(phaseBins)-1
    phaseBinLabel{p} = sprintf('%.3g : %.3g', phaseBins(p), phaseBins(p+1));
end

corrTrlPos = nan(size(trialInfo,1),1);
corrTrlOdr = nan(size(trialInfo,1),1);
for trl = 2:size(trialInfo,1)
    if trialInfo(trl,1)==1 && (trialInfo(trl,3) - trialInfo(trl-1,3) == 1)
        corrTrlPos(trl) = trialInfo(trl,3);
        corrTrlOdr(trl) = trialInfo(trl-1,4);
    end
end

corrTrlEnsmbl = unitEpoch(:,:,~isnan(corrTrlPos));
corrTrlLFP = lfpEpoch(:,:,~isnan(corrTrlPos));
corrTrlPos = corrTrlPos(~isnan(corrTrlPos));
corrTrlOdr = corrTrlOdr(~isnan(corrTrlOdr));     

unitPosFvals = repmat({nan(length(eventTimeBins)-dataBinSize, length(phaseBins)-1, length(bands))}, [1, length(unitIDs)]);
unitOdrFvals = repmat({nan(length(eventTimeBins)-dataBinSize, length(phaseBins)-1, length(bands))}, [1, length(unitIDs)]);

for u = 1:length(unitIDs)
    curTetParts = strsplit(unitIDs{u}, '-');
    tic
    fprintf('Running %s...', unitIDs{u});
    curUnitPosFvals = nan(length(eventTimeBins)-dataBinSize, length(phaseBins)-1, length(bands));
    curUnitOdrFvals = nan(length(eventTimeBins)-dataBinSize, length(phaseBins)-1, length(bands));
    for b = 1:length(bands)
        curLFPcol = strcmp(bands{b}, lfpIDparts(2,:)) & strcmp(curTetParts{1}, lfpIDparts(1,:));
        for p = 1:length(phaseBins)-1
            for t = 1+dataBinSize/2:length(eventTimeBins)-(dataBinSize/2)
                timeMask = false(length(eventTimeBins),1);
                timeMask(t-(dataBinSize/2):t+dataBinSize/2) = true;
                
                spikeLog = corrTrlEnsmbl(timeMask,u,:);
                phaseVals = corrTrlLFP(timeMask,curLFPcol,:);
                phaseValsLog = phaseVals>=phaseBins(p) & phaseVals<phaseBins(p+1);
                
                [~,tablePos,~] = anova1(reshape(sum(spikeLog & phaseValsLog), [size(corrTrlEnsmbl,3),1]), corrTrlPos, 'off');
                if ~isempty(tablePos{2,5})
                    curUnitPosFvals(t,p,b) = tablePos{2,5};
                else
                    curUnitPosFvals(t,p,b) = 1;
                end
                
                [~,tableOdr,~] = anova1(reshape(sum(spikeLog & phaseValsLog), [size(corrTrlEnsmbl,3),1]), corrTrlOdr, 'off');
                if ~isempty(tableOdr{2,5})
                    curUnitOdrFvals(t,p,b) = tableOdr{2,5};
                else
                    curUnitOdrFvals(t,p,b) = 1;
                end
            end
        end
    end
    curUnitPosFvals(1:dataBinSize/2,:,:) = [];
    curUnitPosFvals(isnan(curUnitPosFvals)) = 1;
    unitPosFvals{u} = curUnitPosFvals;
    
    curUnitOdrFvals(1:dataBinSize/2,:,:) = [];
    curUnitOdrFvals(isnan(curUnitOdrFvals)) = 1;
    unitOdrFvals{u} = curUnitOdrFvals;
    toc
end

% Plot this shit.
% subplots = band X 1; plot each phase as a different line plot on the axis
f = figure;
for b = 1:length(bands)
    subplot(length(bands),1,b)
    hold on
    means = nan(1,length(phaseBins)-1);
    for p = 1:length(phaseBins)-1
        curOdrVals = cell2mat(cellfun(@(a)a(:,p,b), unitOdrFvals, 'uniformoutput', 0));
        curPosVals = cell2mat(cellfun(@(a)a(:,p,b), unitPosFvals, 'uniformoutput', 0));
        
        means(p) = plot(timePeriod, mean(curPosVals - curOdrVals,2), 'color', phaseColors{p});
        plot(timePeriod, ((std(curPosVals - curOdrVals,1,2)./(length(unitIDs)-1))+mean(curPosVals - curOdrVals,2)), 'linestyle', ':', 'color', phaseColors{p});
        plot(timePeriod, mean(curPosVals - curOdrVals,2)-(std(curPosVals - curOdrVals,1,2)./(length(unitIDs)-1)), 'linestyle', ':', 'color', phaseColors{p});
    end
    title(sprintf('Information Bias: Current Position vs Previous Odor (%s) (%ims window)', bands{b},dataBinSize), 'interpreter', 'none')
    line(get(gca, 'xlim'), [0 0], 'color', 'k', 'linewidth', 1);
    set(gca, 'ylim', [-0.5 1]);
    grid on;
end
legend(means, phaseBinLabel);
annotation(f,'textbox', [0.01 0.01 0.96 0.03], 'FitBoxToText','off', 'string', cd, 'interpreter', 'none', 'linestyle', 'none');
orient(f, 'tall');
% print(f);

% Plot this shit per unit
for u = 1:length(unitIDs)
    f = figure;
    sps = 1:length(bands);
    for b = 1:length(bands)
        sps(b) = subplot(length(bands),1,b);
        hold on
        means = nan(1,length(phaseBins)-1);
        for p = 1:length(phaseBins)-1
            curOdrVals = unitOdrFvals{u}(:,p,b);
            curPosVals = unitPosFvals{u}(:,p,b);
            
            means(p) = plot(timePeriod', curPosVals - curOdrVals);%, 'color', phaseColors{p});
        end
        title(sprintf('Information Bias: Current Position vs Previous Odor (%s) (%ims window)', bands{b},dataBinSize), 'interpreter', 'none')
%         line(get(gca, 'xlim'), [0 0], 'color', 'k', 'linewidth', 1);
%         set(gca, 'ylim', [-0.5 1]);
        grid on;
    end
    linkaxes(sps, 'y');
    fnlSpsPos = get(sps(end), 'position');
    h = legend(means, phaseBinLabel);
    set(h, 'position', [fnlSpsPos(1)+fnlSpsPos(3)+0.01 fnlSpsPos(2) 0.05 fnlSpsPos(4)]);
    annotation(f,'textbox', [0.01 0.01 0.96 0.03], 'FitBoxToText','off', 'string', sprintf('%s: %s', unitIDs{u}, cd), 'interpreter', 'none', 'linestyle', 'none');
    drawnow
    orient(f, 'tall');
%     print(f)
end
