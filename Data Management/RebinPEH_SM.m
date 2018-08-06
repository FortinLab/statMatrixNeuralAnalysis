function [meanPEH, newBins] = RebinPEH_SM(trialOrgSpikes, window, newBinSize)

if iscell(trialOrgSpikes)
    spkLog = logical(cell2mat(trialOrgSpikes)');
elseif islogical(trialOrgSpikes)
    spkLog = trialOrgSpikes';
else
    spkLog = logical(trialOrgSpikes)';
end

oldBins = linspace(window(1), window(2), size(spkLog,2));
oldBinMtx = repmat(oldBins, [size(spkLog,1),1]);

oldSpkRelTimes = oldBinMtx(spkLog);

% newBins = linspace(window(1), window(2), (sum(abs(window))/newBinSize)+1);
newBins = window(1):newBinSize:window(2);
newHist = histcounts(oldSpkRelTimes, newBins);

meanPEH = newHist./size(spkLog,1);



