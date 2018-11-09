

origDir = cd;
[fileDir] = uigetdir(origDir);
if fileDir==0
    disp('Analysis Cancelled')
    return
else
    cd(fileDir)
end


%% Define Parameters
alignment = 'PokeOut';
windowStart = -1;
windowEnd = 1;
dataBinSize = 125;
%% Load Relevant Data
% 'trialInfo' : Matrix containing information about the identity and
%           outcome of each trial. Rows represent individual trials;
%           columns are organized thusly:
%               Column 1) Trial performance. 1 = Correct, 0 = Incorrect
%               Column 2) Trial InSeq log. 1 = InSeq, 0 = OutSeq
%               Column 3) Trial Position.
%               Column 4) Trial Odor.
%               Column 5) Poke Duration.
[unitEpoch, unitIDs, ~, ~, ~, eventTimeBins, trialInfo] = EpochExtraction_SM(alignment, windowStart, windowEnd, 'org', 'TiUTr', 'lfpBand', 'Beta');

%% Bin the unit data
samp = mode(diff(eventTimeBins));
binLims = linspace(windowStart, windowEnd, (windowEnd-windowStart)/(dataBinSize*samp)+2);
binnedUnitEpoch = nan(length(binLims)-1, size(unitEpoch,2), size(unitEpoch,3));
for uni = 1:length(unitIDs)
    for trl = 1:size(trialInfo,1)
        binnedUnitEpoch(:,uni,trl) = histcounts(eventTimeBins(unitEpoch(:,uni,trl)==1), binLims);
    end
end
binnedTimeBins = binLims(2:end)-(dataBinSize*samp/2);

%%
rt = [ssnData.ReactionTime];
rt(isnan(rt))= 0;
perfLog = trialInfo(:,1)==1;
isLog = trialInfo(:,2)==1;
aLog = (trialInfo(:,3)==1 | trialInfo(:,4)==1);
% Select data
corrTrlEnsmbl= binnedUnitEpoch(:,:,isLog & perfLog);

corrVals = nan(length(binnedTimeBins), length(unitIDs));
sigVals = nan(length(binnedTimeBins), length(unitIDs));
for u = 1:length(unitIDs)
    for t = 1:length(binnedTimeBins) 
        [r,p] = corrcoef(reshape(corrTrlEnsmbl(t,u,:), [size(corrTrlEnsmbl,3),1]), rt(isLog & perfLog));
        corrVals(t,u) = r(2);
        sigVals(t,u) = p(2);
    end
end
corrVals(isnan(corrVals)) = 0;

figure
plot(binnedTimeBins, corrVals)

%%
for u = 1:length(unitIDs)
    figure
    a = subplot(2,1,1);
    plot(binnedTimeBins, corrVals(:,u), 'color', 'black');
    xlabel(sprintf('Time relative to %s(s)', alignment));
    title([unitIDs{u} ' r-Values']);
    set(gca, 'ylim', [-1 1]);
    b = subplot(2,1,2);
    plot(binnedTimeBins, sigVals(:,u), 'color', 'black');
    xlabel(sprintf('Time relative to %s(s)', alignment));
    title([unitIDs{u} ' p-Values']);
    hold on;
    line([-1 1], [0.05 0.05], 'linestyle', '--', 'color', 'black');
    set(gca, 'ylim', [0 0.5]);
    linkaxes([a,b], 'x');
    drawnow
end