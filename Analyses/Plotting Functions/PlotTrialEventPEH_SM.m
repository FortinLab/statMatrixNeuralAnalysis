function PlotTrialEventPEH_SM(unitID, behavMatrices, behavMatrixIDs, groupingLogs, groupingLogIDs, curUniSpikeLog, origBinWindows, pehBinSize, saveYN)
%% PlotTrialEventPEH_SM
% Future versions should have:
%   1) Clear description up here instead of a to-do list
%   2) Code should be condensed down to a single for loop.
%       - To do this the inputs will need to be changed to be separated
%           into: 
%               : BehaviorMatrices
%               : TrialSelectionVectors (currently referred to as "eventLog1/2")
%               : Spike Vector 
%               : Original Event Windows
%               : Figure ID
%
%% 
figure('Name', sprintf('%s (%s vs %s)', unitID, groupingLogIDs{end-1}, groupingLogIDs{end}), 'NumberTitle', 'off');
annotation('textbox', [0.05 0.9 0.9 0.1], 'String', sprintf('%s (%s vs %s)', unitID, groupingLogIDs{end-1}, groupingLogIDs{end}), 'linestyle', 'none');

%% 
if isstruct(groupingLogs)
    groupingLogs = cell2mat(struct2cell(groupingLogs));
    groupingLogIDs = fieldnames(groupingLogs);
elseif iscell(groupingLogs)
    groupingLogs = cell2mat(groupingLogs');
elseif islogical(groupingLogs)
end


numSubplots = length(behavMatrixIDs) * length(groupingLogIDs);
subplotKey = reshape(1:numSubplots, [length(behavMatrixIDs) length(groupingLogIDs)]);

%% Poke In Aligned Plots
subplotIDs = nan(1,numSubplots);
for eve = 1:length(behavMatrixIDs)
    curEvent = behavMatrices(eve,:);
    curEventID = behavMatrixIDs{eve};
    curEveSubplots = subplotKey(eve,:);
    for grp = 1:length(groupingLogIDs)
        curGroupLog = groupingLogs(grp,:);
        curGroupID = groupingLogIDs{grp};       
        subplotIDs(curEveSubplots(grp)) = subplot(length(groupingLogIDs),length(behavMatrixIDs), curEveSubplots(grp));
        
        curEventData = ExtractTrialData_SM(curEvent(curGroupLog), curUniSpikeLog);
        noSpkLog = cellfun(@(a)isempty(a), curEventData);
        curEventData(noSpkLog) = [];
        curEventPEH = cell(length(curEventData),1);
        if ~isempty(curEventData)
            for trl = 1:length(curEventPEH)
                [curEventPEH{trl,1}, newBins] = RebinPEH_SM(curEventData{trl}, origBinWindows, pehBinSize);
            end
            meanEventPEH = mean(cell2mat(curEventPEH));   
            semEventPEH = std(cell2mat(curEventPEH))./sqrt(length(curEventPEH)-1);
            BarPlotErrorbars(meanEventPEH,semEventPEH, 'Color', 'Black', 'XTick', newBins(1:end-1)+0.05);
            axis tight
        else
            set(subplotIDs(curEveSubplots(grp)), 'xlim', origBinWindows, 'ylim', [0 0.0001]);
        end
        title(sprintf('%s: %s', curEventID, curGroupID));
        drawnow
    end
end
%% Link Axes!
linkaxes(subplotIDs, 'xy');

if saveYN==1
    set(figID, 'PaperOrientation', 'landscape');
    print('-fillpage', figID, '-dpdf', sprintf('%s (%s vs %s).pdf', unitID, event1ID, event2ID));
end