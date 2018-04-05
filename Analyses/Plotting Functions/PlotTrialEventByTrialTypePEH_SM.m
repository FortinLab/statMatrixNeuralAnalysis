function PlotTrialEventByTrialTypePEH_SM(unitID, behavMatrices, behavMatrixIDs, trialTypeLog, trialTypeIDs, groupingLogs, groupingLogIDs, curUniSpikeLog, origBinWindows, pehBinSize)
%% PlotTrialEventByTrialTypePEH_SM
%   Creates new figures for each trial type
%   Should be updated to be more flexible and have subplots dynamically
%   assigned, ideally with a single for loop. For now, hard coded all the
%   way.

%%
if isstruct(trialTypeLog)
    trialTypes = cell2mat(struct2cell(trialTypeLog));
    trialTypeIDs = fieldnames(trialTypeLog);
elseif iscell(trialTypeLog)
    trialTypes = cell2mat(trialTypeLog');
elseif islogical(trialTypeLog)
    trialTypes = trialTypeLog;
end

if isstruct(groupingLogs)
    groupingLogs = cell2mat(struct2cell(groupingLogs));
    groupingLogIDs = fieldnames(groupingLogs);
elseif iscell(groupingLogs)
    groupingLogs = cell2mat(groupingLogs');
elseif islogical(groupingLogs)
end

numSubplots = length(trialTypeIDs) * length(groupingLogIDs);
subplotKey = reshape(1:numSubplots, [length(groupingLogIDs) length(trialTypeIDs)]);

for eve = 1:size(behavMatrices,1)
    curEvent = behavMatrices(eve,:);
    curEventID = behavMatrixIDs{eve};
    figure('Name', sprintf('%s %s (%s vs %s)', unitID, curEventID, groupingLogIDs{end-1}, groupingLogIDs{end}), 'NumberTitle', 'off');
    annotation('textbox', [0.05 0.9 0.9 0.1], 'String', sprintf('%s %s (%s vs %s)', unitID, curEventID, groupingLogIDs{end-1}, groupingLogIDs{end}), 'linestyle', 'none');
    subplotIDs = nan(1,numSubplots);
    for tt = 1:length(trialTypeIDs)
        curTTLog = trialTypes(tt,:);
        curTTid = trialTypeIDs{tt};
        curTTsubplots = subplotKey(:,tt);
        for grp = 1:size(groupingLogs,1)
            curGroupLog = groupingLogs(grp,:);
            curGroupID = groupingLogIDs{grp};
            subplotIDs(curTTsubplots(grp)) = subplot(length(trialTypeIDs),length(groupingLogIDs), curTTsubplots(grp));
                        
            curEventData = ExtractTrialData_SM(curEvent(curTTLog & curGroupLog), curUniSpikeLog);
            noSpkLog = cellfun(@(a)isempty(a), curEventData);
            curEventData(noSpkLog) = [];
            if ~isempty(curEventData)
                [curEventPEH, newBins] = RebinPEH_SM(curEventData, origBinWindows, pehBinSize);
                bar(newBins(1:end-1)+0.05, curEventPEH);
                axis tight
            else
                set(subplotIDs(curTTsubplots(grp)), 'xlim', [0 0.00001], 'ylim', [0 0.000001]);
            end
            title(sprintf('%s %s: %s', curTTid, curEventID, curGroupID));
        end
    end
    linkaxes(subplotIDs, 'xy');
    ylims = get(subplotIDs(end), 'ylim');
    if ylims(1)<0
        set(subplotIDs(end), 'ylim', [0 0.25]);
    end
end

