function [binnedData, newBins] = RebinPEH2_SM(trialOrgSpikes, window, newBinSize, sampleRate)


if iscell(trialOrgSpikes)
    spkLog = logical(cell2mat(trialOrgSpikes)');
elseif islogical(trialOrgSpikes)
    spkLog = trialOrgSpikes';
else
    spkLog = logical(trialOrgSpikes)';
end

if ~isempty(window)
    newBins = window(1):newBinSize:window(2);
    oldBins = linspace(window(1), window(2), size(spkLog,2));
else
    if nargin==3 || isempty(sampleRate)
        error('If not specifying a time period for the window, you must specify a samplerate at which the data is collected/binned');
    else
        windowSize = (size(spkLog,2)*1)/sampleRate;
        if windowSize<newBinSize
            newBinSize = windowSize-(1/sampleRate);
        end
        newBins = 0:newBinSize:windowSize;
        oldBins = (0:size(spkLog,2)-1)/sampleRate;
    end
end
binnedData = nan(size(spkLog,1), length(newBins)-1);
for row = 1:size(spkLog,1)
    oldSpkRelTimes = oldBins(spkLog(row,:));

    % newBins = linspace(window(1), window(2), (sum(abs(window))/newBinSize)+1);
    binnedData(row,:) = histcounts(oldSpkRelTimes, newBins);
end



