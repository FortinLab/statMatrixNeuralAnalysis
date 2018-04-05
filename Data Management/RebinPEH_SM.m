function [meanPEH, newBins] = RebinPEH_SM(trialOrgSpikes, window, newBinSize)

spkLog = logical(cell2mat(trialOrgSpikes)');

oldBins = linspace(window(1), window(2), length(trialOrgSpikes{1}));
oldBinMtx = repmat(oldBins, [size(spkLog,1),1]);

oldSpkRelTimes = oldBinMtx(spkLog);

newBins = linspace(window(1), window(2), (sum(abs(window))/newBinSize)+1);
newHist = histcounts(oldSpkRelTimes, newBins);

meanPEH = newHist./length(trialOrgSpikes);



