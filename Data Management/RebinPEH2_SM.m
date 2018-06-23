function [binnedData, newBins] = RebinPEH2_SM(trialOrgSpikes, window, newBinSize)

if iscell(trialOrgSpikes)
    spkLog = logical(cell2mat(trialOrgSpikes)');
elseif islogical(trialOrgSpikes)
    spkLog = trialOrgSpikes';
else
    spkLog = logical(trialOrgSpikes)';
end

newBins = window(1):newBinSize:window(2);
oldBins = linspace(window(1), window(2), size(spkLog,2));
binnedData = nan(size(spkLog,1), length(newBins)-1);
for row = 1:size(spkLog,1)
    oldSpkRelTimes = oldBins(spkLog(row,:));

    % newBins = linspace(window(1), window(2), (sum(abs(window))/newBinSize)+1);
    binnedData(row,:) = histcounts(oldSpkRelTimes, newBins);
end



