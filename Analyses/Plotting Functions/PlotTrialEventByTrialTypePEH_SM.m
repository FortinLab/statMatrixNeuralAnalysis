function [subplotData] = PlotTrialEventByTrialTypePEH_SM(unitID, behavMatrices, behavMatrixIDs, trialTypeLog, trialTypeIDs, groupingLogs, groupingLogIDs, curUniSpikeLog, origBinWindows, pehBinSize, saveYN)
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
[ttIDstart, ttIDend] = regexp(trialTypeIDs{1}, '([A-Za-z]*)');
ttID = trialTypeIDs{1}(ttIDstart:ttIDend);

if isstruct(groupingLogs)
    groupingLogs = cell2mat(struct2cell(groupingLogs));
    groupingLogIDs = fieldnames(groupingLogs);
elseif iscell(groupingLogs)
    groupingLogs = cell2mat(groupingLogs');
elseif islogical(groupingLogs)
end

numSubplots = length(trialTypeIDs) * length(groupingLogIDs);
subplotKey = reshape(1:numSubplots, [length(groupingLogIDs) length(trialTypeIDs)]);
figIDs = nan(1,size(behavMatrices,1));

subplotIDs = nan(size(behavMatrices,1),numSubplots);
subplotData = [];
for eve = 1:size(behavMatrices,1)
    curEvent = behavMatrices(eve,:);
    curEventID = behavMatrixIDs{eve};
    figIDs(eve) = figure('Name', sprintf('%s %s (%s vs %s)', unitID, curEventID, groupingLogIDs{end-1}, groupingLogIDs{end}), 'NumberTitle', 'off');
    annotation('textbox', [0.05 0.9 0.9 0.1], 'String', sprintf('%s %s (%s vs %s)', unitID, curEventID, groupingLogIDs{end-1}, groupingLogIDs{end}), 'linestyle', 'none');
    for tt = 1:length(trialTypeIDs)
        curTTLog = trialTypes(tt,:);
        curTTid = trialTypeIDs{tt};
        curTTsubplots = subplotKey(:,tt);
        for grp = 1:size(groupingLogs,1)
            curGroupLog = groupingLogs(grp,:);
            curGroupID = groupingLogIDs{grp};
            subplotIDs(eve,curTTsubplots(grp)) = subplot(length(trialTypeIDs),length(groupingLogIDs), curTTsubplots(grp));
                        
            curEventData = ExtractTrialData_SM(curEvent(curTTLog & curGroupLog), curUniSpikeLog);
            noSpkLog = cellfun(@(a)isempty(a), curEventData);
            curEventData(noSpkLog) = [];
            curEventPEH = cell(length(curEventData),1);
            if ~isempty(curEventData)
                for trl = 1:length(curEventPEH)
                    [curEventPEH{trl,1}, newBins] = RebinPEH_SM(curEventData{trl}, origBinWindows, pehBinSize);
                end
                meanEventPEH = mean(cell2mat(curEventPEH),1);
                semEventPEH = std(cell2mat(curEventPEH),0,1)./sqrt(length(curEventPEH)-1);
                BarPlotErrorbars(meanEventPEH, semEventPEH, 'Color', 'Black', 'XTick', newBins(1:end-1)+(mode(diff(newBins))/2));
                axis tight
            else
                set(subplotIDs(eve,curTTsubplots(grp)), 'xlim', origBinWindows, 'ylim', [0 0.0001]);
            end
            title(sprintf('%s %s: %s Trials', curTTid, curEventID, curGroupID));
            subplotData.(curEventID).(curGroupID).(curTTid) = curEventData;
        end
    end
    linkaxes(subplotIDs(eve,:), 'xy');
    ylims = get(subplotIDs(eve,end), 'ylim');
    if ylims(1)<0
        set(subplotIDs(eve,end), 'ylim', [0 ylims(2)]);
    end    
end
linkaxes(subplotIDs, 'xy');

if saveYN==1
    for fig = 1:length(figIDs)
        set(figIDs(fig), 'PaperOrientation', 'landscape');
        print('-fillpage', figIDs(fig), '-dpdf', sprintf('%s %s by %s (%s vs %s)', unitID, behavMatrixIDs{fig}, ttID, groupingLogIDs{end-1}, groupingLogIDs{end}));
    end
end



