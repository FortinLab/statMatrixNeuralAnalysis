function PlotSpikePhaseRelations_SM(statMatrix, statMatrixColIDs, behavMatrices, trialTypeLog, trialTypeIDs, groupingLogs, groupingLogIDs, trialWindows, saveYN)
%% Identify Input Variables
% Determine Trial Types
if isstruct(trialTypeLog)
    trialTypes = cell2mat(struct2cell(trialTypeLog));
    trialTypeIDs = fieldnames(trialTypeLog);
elseif iscell(trialTypeLog)
    trialTypes = cell2mat(trialTypeLog');
elseif islogical(trialTypeLog)
    trialTypes = trialTypeLog;
end
[ttIDstart, ttIDend] = regexp(trialTypeIDs{1}, '([A-Za-z]*)');
ttID = trialTypeIDs{1}(ttIDstart:ttIDend);

% Determine the grouping factors
if isstruct(groupingLogs)
    groupingLogs = cell2mat(struct2cell(groupingLogs));
    groupingLogIDs = fieldnames(groupingLogs);
elseif iscell(groupingLogs)
    groupingLogs = cell2mat(groupingLogs');
elseif islogical(groupingLogs)
end

% Identify the LFP bands filtered
phaseValsLog = cellfun(@(a)~isempty(a), regexp(statMatrixColIDs, '_HilbVals\>'));
rawChanLog = cellfun(@(a)~isempty(a), strfind(statMatrixColIDs, 'Raw'));
phaseValCols = find(phaseValsLog & ~rawChanLog);
lfpBands = cellfun(@(a)a{3}, cellfun(@(a)strsplit(a, '_'), statMatrixColIDs(phaseValCols), 'uniformoutput', 0), 'uniformoutput', 0);

% Identify units
uniLog = cellfun(@(a)~isempty(a), strfind(statMatrixColIDs, '-U'));
unitIDs = statMatrixColIDs(uniLog);
unitCols = find(uniLog);

%% Define Plotting Constants
phaseBins = -pi:pi/8:pi;
timeBins = trialWindows(1):0.1:trialWindows(2);

%% Plot Trial Data by Band
phaseOverTime = repmat({nan(length(trialTypeIDs), length(groupingLogIDs),length(lfpBands))}, [2,sum(uniLog)]);
for uni = 1:sum(uniLog)
    for bnd = 1:length(lfpBands)
        figure;
        colPlots = nan([length(groupingLogIDs), length(trialTypeIDs)]);
        cLimVals = cell(length(lfpBands));
        annotation('textbox', [0.05 0.9 0.9 0.1], 'String', sprintf('%s %s Phase Relations (%s vs %s)', unitIDs{uni}, lfpBands{bnd}, groupingLogIDs{1}, groupingLogIDs{2}), 'linestyle', 'none');
        for tt = 1:length(trialTypeIDs)
            for grp = 1:length(groupingLogIDs)
                curTrls = behavMatrices(trialTypes(tt,:) & groupingLogs(grp,:));
                tempPhaseVals = cell(length(curTrls),1);
                tempTimeVals = cell(length(curTrls),1);
                for trl = 1:length(curTrls)
                    trlTSs = statMatrix(curTrls(trl).TrialLogVect,1);
                    tempPhaseVals{trl} = statMatrix(curTrls(trl).TrialLogVect & statMatrix(:,unitCols(uni))>=1, phaseValCols(bnd));
                    tempTimeVals{trl} = statMatrix(curTrls(trl).TrialLogVect & statMatrix(:,unitCols(uni))>=1, 1) - trlTSs(1) + trialWindows(1);
                end
                sp = sub2ind([length(groupingLogIDs), length(trialTypeIDs)], grp, tt);
                colPlots(sp) = subplot(length(trialTypeIDs), length(groupingLogIDs), sp);
%                 scatter(cell2mat(tempTimeVals), cell2mat(tempPhaseVals), '*k');
%                 set(gca, 'xlim', trialWindows, 'ylim', [-pi pi]);
                contourf(timeBins(1:end-1), phaseBins(1:end-1),(histcounts2(cell2mat(tempTimeVals), cell2mat(tempPhaseVals), timeBins, phaseBins)')./length(curTrls),20, 'linestyle', 'none');
%                 imagesc(timeBins, phaseBins, flipud(histcounts2(cell2mat(tempTimeVals), cell2mat(tempPhaseVals), timeBins, phaseBins)'));
                cLimVals{sp} = get(colPlots(sp), 'CLim');
                title([trialTypeIDs{tt} ' ' groupingLogIDs{grp}], 'interpreter', 'none');
%                 drawnow;
            end
        end
        cLimMax = max(median(cell2mat(cLimVals(:))));
        for sp = 1:length(colPlots(:))
            set(colPlots(sp), 'clim', [0 cLimMax]);
        end
        colormap jet
        if saveYN==1
            set(gcf, 'PaperOrientation', 'landscape');
            print('-fillpage', gcf, '-dpdf', sprintf('%s %s Phase Relations (%s vs %s)', unitIDs{uni}, lfpBands{bnd}, groupingLogIDs{1}, groupingLogIDs{2}));
        end
    end
end

