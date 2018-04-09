function PlotTrialEventSpect_SM(unitID, behavMatrices, behavMatrixIDs, groupingLogs, groupingLogIDs, statMatrix, statMatrixColIDs, eventWindow, freqWindow, saveYN)
%% PlotTrialEventSpect_SM
%   Put something clever and descriptive here...
%
%%
figure('Name', sprintf('%s (%s vs %s)', unitID, groupingLogIDs{end-1}, groupingLogIDs{end}), 'NumberTitle', 'off');
annotation('textbox', [0.05 0.9 0.9 0.1], 'String', sprintf('%s (%s vs %s)', unitID, groupingLogIDs{end-1}, groupingLogIDs{end}), 'linestyle', 'none');

%% Pull out Raw LFP
rawLFPcolLog = cellfun(@(a)~isempty(a), regexp(statMatrixColIDs, 'LFP_Raw$'));
rawLFP = statMatrix(:,rawLFPcolLog);

sampleRate = 1/mode(diff(statMatrix(:,1)));

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

%%
subplotIDs = nan(1,numSubplots);
subplotClims = cell(numSubplots,1);
for eve = 1:length(behavMatrixIDs)
    curEvent = behavMatrices(eve,:);
    curEventID = behavMatrixIDs{eve};
    curEveSubplots = subplotKey(eve,:);
    for grp = 1:length(groupingLogIDs)
        curGroupLog = groupingLogs(grp,:);
        curGroupID = groupingLogIDs{grp};
       
        subplotIDs(curEveSubplots(grp)) = subplot(length(groupingLogIDs),length(behavMatrixIDs), curEveSubplots(grp));
        curEventLFP = ExtractTrialData_SM(curEvent(curGroupLog), rawLFP);
        noLFPlog = cellfun(@(a)isempty(a), curEventLFP);
        curEventLFP(noLFPlog) = [];
        curEventSpect = cell(1,1,length(curEventLFP));
        if ~isempty(curEventLFP)
            for trl = 1:length(curEventLFP)
                curEventSpect{1,1,trl} = MorletAG(curEventLFP{trl}, sampleRate, freqWindow(1), freqWindow(2));
            end
%             meanEventSpect = mean(cell2mat(curEventSpect),3);            
            meanEventSpect = median(cell2mat(curEventSpect),3);            
            contourf(eventWindow(1):1/sampleRate:eventWindow(2), freqWindow(1):freqWindow(2), meanEventSpect, 20, 'linestyle', 'none');
            subplotClims{curEveSubplots(grp)} = get(subplotIDs(curEveSubplots(grp)), 'clim');
        else
            set(subplotIDs(curEveSubplots(grp)), 'xlim', eventWindow, 'ylim', freqWindow);
        end
        title(sprintf('%s: %s', curEventID, curGroupID));
        drawnow
    end
end
allClims = cell2mat(subplotClims);
commonClims = [median(allClims(:,1)) median(allClims(:,2))];
for sp = 1:numSubplots
    set(subplotIDs(sp), 'clim', commonClims);
end

if saveYN==1
    set(gcf, 'PaperOrientation', 'landscape');
    print('-fillpage', gcf, '-dpdf', sprintf('%s Spectrogram (%s vs %s)', unitID, groupingLogIDs{end-1}, groupingLogIDs{end}));
end

