function PositionInformationContent_SM

origDir = cd;
[fileDir] = uigetdir(origDir);
if fileDir==0
    disp('Analysis Cancelled')
    return
else
    cd(fileDir)
end
dirContents = dir(fileDir);
fileNames = {dirContents.name};

%% Define Standard Variables
timeBinSize = 0.025;
slideWindowSize = 30;
spatialBinSize = 2;

%% Load Relevant Data
load(fileNames{cellfun(@(a)~isempty(a), strfind(fileNames, 'BehaviorMatrix'))});
load(fileNames{cellfun(@(a)~isempty(a), strfind(fileNames, 'EnsembleMatrix'))});
load(fileNames{cellfun(@(a)~isempty(a), strfind(fileNames, 'OrientationMatrix'))});

%% Pre-Process the Orientation Matrix Data
portPosX = strcmp(orientMatrixColIDs, 'PortX');
portPosY = strcmp(orientMatrixColIDs, 'PortY');
headPosX = strcmp(orientMatrixColIDs, 'HeadX');
headPosY = strcmp(orientMatrixColIDs, 'HeadY');
tailPosX = strcmp(orientMatrixColIDs, 'TailX');
tailPosY = strcmp(orientMatrixColIDs, 'TailY');

portAngle = nan(size(orientMatrix,1),1);
headAngle = nan(size(orientMatrix,1),1);
tailAngle = nan(size(orientMatrix,1),1);

htVal = nan(size(orientMatrix,1),1);
hpVal = nan(size(orientMatrix,1),1);
ptVal = nan(size(orientMatrix,1),1);

posIndices = find(~isnan(orientMatrix(:,2)));

for pos = 1:length(posIndices)
    curPosNdx = posIndices(pos);
    curPortX = orientMatrix(curPosNdx, portPosX);
    curPortY = orientMatrix(curPosNdx, portPosY);
    curHeadX = orientMatrix(curPosNdx, headPosX);
    curHeadY = orientMatrix(curPosNdx, headPosY);
    curTailX = orientMatrix(curPosNdx, tailPosX);
    curTailY = orientMatrix(curPosNdx, tailPosY);
    
    htVal(curPosNdx) = sqrt((curTailX - curHeadX)^2 + (curTailY - curHeadY)^2);
    hpVal(curPosNdx) = sqrt((curPortX - curHeadX)^2 + (curPortY - curHeadY)^2);
    ptVal(curPosNdx) = sqrt((curTailX - curPortX)^2 + (curTailY - curPortY)^2);
    
    portAngle(curPosNdx) = rad2deg(acos((ptVal(curPosNdx)^2 + hpVal(curPosNdx)^2 - htVal(curPosNdx)^2)/(2*ptVal(curPosNdx)*hpVal(curPosNdx))));
    headAngle(curPosNdx) = rad2deg(acos((htVal(curPosNdx)^2 + hpVal(curPosNdx)^2 - ptVal(curPosNdx)^2)/(2*htVal(curPosNdx)*hpVal(curPosNdx))));
    tailAngle(curPosNdx) = rad2deg(acos((htVal(curPosNdx)^2 + ptVal(curPosNdx)^2 - hpVal(curPosNdx)^2)/(2*htVal(curPosNdx)*ptVal(curPosNdx))));    
end

orientMatrix = [orientMatrix, portAngle, headAngle, tailAngle, htVal, hpVal, ptVal];
orientMatrixColIDs = [orientMatrixColIDs, {'PortAngle'}, {'HeadAngle'}, {'TailAngle'}, {'HeadTailLength'}, {'HeadPortLength'}, {'PortTailLength'}];

%% Calculate Instantaneous Firing Rates
% For this I'm opting to use the instantaneous firing rate (i.e. convolved
% with a gaussian) because the camera capture rate is not exactly uniform,
% which is what is needed for using discrete time bins. This may lead to
% some distortion in the spiking value... I may want to do it with discrete
% bins down the road to verify I get similar values.

% Create Gaussian
instFRgauss = gausswin(slideWindowSize);
instFRgauss = instFRgauss/(length(instFRgauss)*mode(diff(behavMatrix(:,1))));

% Convolve Gaussian with Binary Spike Trains
uniInstFR = [ensembleMatrix(:,1), nan(size(ensembleMatrix,1), size(ensembleMatrix,2)-1)];
for uni = 2:size(ensembleMatrix,2)
    uniInstFR(:,uni) = conv(ensembleMatrix(:,uni), instFRgauss, 'same');
end

%% Organize Trial Data
% Create Extraction Matrices
pokeInTrialMatrix = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [-0.25 0.45], 'PokeIn');
pokeOutTrialMatrix = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [-0.45 0.25], 'PokeOut');

% Create Trial Logical Vectors
corrTrlLog = [pokeInTrialMatrix.Performance];
isLog = [pokeInTrialMatrix.TranspositionDistance]==0;
pokeInNdxs = {pokeInTrialMatrix.PokeInIndex};
pokeInNdxsISC = pokeInNdxs(corrTrlLog & isLog);
odorIDsISC = [pokeInTrialMatrix(corrTrlLog & isLog).Odor];

% Extract Orientation Data
% Pull out the InSeq Correct Trials
pokeInOrientationAll = ExtractTrialData_SM(pokeInTrialMatrix, orientMatrix); %#ok<*NODEF>
pokeInOrientationISC = pokeInOrientationAll(corrTrlLog & isLog);
pokeOutOrientationAll = ExtractTrialData_SM(pokeOutTrialMatrix, orientMatrix);
pokeOutOrientationISC = pokeOutOrientationAll(corrTrlLog & isLog);
% Now concatenate them into a single cell vector
orientationTrialDataRaw = cellfun(@(a,b)[a;b], pokeInOrientationISC, pokeOutOrientationISC, 'uniformoutput', 0);
% Now correct the timestamps in the vector so that poke In and poke Out are
% equi-distant.
relTimeVects = cellfun(@(a,b,c)[a(:,1); b(:,1) - (b(1,1)-a(end,1)) + (b(2,1)-b(1,1))] - orientMatrix(c,1), pokeInOrientationISC, pokeOutOrientationISC, pokeInNdxsISC, 'uniformoutput',0);
orientationTrialData = cellfun(@(a,b)[a,b(:,2:end)], relTimeVects, orientationTrialDataRaw, 'uniformoutput', 0);

% Extract Spike Data
pokeInFRall = ExtractTrialData_SM(pokeInTrialMatrix, uniInstFR); %#ok<*NODEF>
pokeInFRisc = pokeInFRall(corrTrlLog & isLog);
pokeOutFRall = ExtractTrialData_SM(pokeOutTrialMatrix, uniInstFR);
pokeOutFRisc = pokeOutFRall(corrTrlLog & isLog);
% Now concatenate them into a single cell vector
firingRateTrialDataRaw = cellfun(@(a,b)[a;b], pokeInFRisc, pokeOutFRisc, 'uniformoutput', 0);
relTimeVects = cellfun(@(a,b,c)[a(:,1); b(:,1) - (b(1,1)-a(end,1)) + (b(2,1)-b(1,1))] - orientMatrix(c,1), pokeInFRisc, pokeOutFRisc, pokeInNdxsISC, 'uniformoutput',0);
firingRateTrialData = cellfun(@(a,b)[a,b(:,2:end)], relTimeVects, firingRateTrialDataRaw, 'uniformoutput', 0);

% Create Trial Position cell vectors
trialIDall = arrayfun(@(a,b)[ones(sum(a.TrialLogVect),1)*a.Odor; ones(sum(b.TrialLogVect),1)*b.Odor], pokeInTrialMatrix, pokeOutTrialMatrix, 'uniformoutput', 0);
trialIDvals = trialIDall(corrTrlLog & isLog);

%% Calculate Standard Variables Used in Analysis
% Average Firing Rate
% For this I need to extract the firing rate values for each unit every
% time there's a position value taken. This is simply because I'm using the
% instantaneous firing rate for expedience. If I were simply binning the
% data into discrete time bins then there wouldn't be an issue, I would
% just calculate the average across bins. Since I'm not, I need to extract
% the instantaneous rate at each capture instance and compute the average
% from that.
% I'm using cellfun here to save space/time.
%   I'm using the trial wise orientation data to create a logical vector of
%   capture times and then applying it to the spike date. Since they're
%   both organized by trial with the same indices extracted it's simple to
%   do.
allFRvals = cell2mat(cellfun(@(a,b)b(~isnan(a(:,2)),2:end), orientationTrialData, firingRateTrialData, 'uniformoutput', 0)');
allTimeVals = cell2mat(cellfun(@(a,b)b(~isnan(a(:,2)),1), orientationTrialData, firingRateTrialData, 'uniformoutput', 0)');
meanFRvals = mean(allFRvals,1);

% Pull out the trial position value to match each of the camera capture
% times
allTrialIDvals = cell2mat(cellfun(@(a,b)b(~isnan(a(:,2)),1), orientationTrialData, trialIDvals, 'uniformoutput', 0)');

% Orientation Occupancy Bins: Define spaces and calculate p(occupy)
% First break the orientation data out of the trial organization.
% Note: 1pixel = ~1.7mm
allOrientData = cell2mat(cellfun(@(a)a(~isnan(a(:,2)),:), orientationTrialData, 'uniformoutput', 0)');
% Head Position: Head X/Y
headXcolLog = strcmp(orientMatrixColIDs, 'HeadX');
headXbins = 30:spatialBinSize:90;
headYcolLog = strcmp(orientMatrixColIDs, 'HeadY');
headYbins = 240:spatialBinSize:300;
headOccupancyMatrix = histcounts2(allOrientData(:,headXcolLog), allOrientData(:,headYcolLog), headXbins, headYbins, 'Normalization', 'probability')';
headOccupancyMatrix(headOccupancyMatrix==0) = nan;
% Tail Position: Tail X/Y
tailXcolLog = strcmp(orientMatrixColIDs, 'TailX');
tailXbins = 60:spatialBinSize:180;
tailYcolLog = strcmp(orientMatrixColIDs, 'TailY');
tailYbins = 180:spatialBinSize:264;
tailOccupancyMatrix = histcounts2(allOrientData(:,tailXcolLog), allOrientData(:,tailYcolLog), tailXbins, tailYbins, 'Normalization', 'probability')';
tailOccupancyMatrix(tailOccupancyMatrix==0) = nan;
% Head Angle: Angle at the head between the port and the tail
headAngleColLog = strcmp(orientMatrixColIDs, 'HeadAngle');
headAngleBins = 0:15:180;
headAngleOccupancyVector = histcounts(allOrientData(:,headAngleColLog), headAngleBins, 'Normalization', 'probability');
% Tail Angle: Angle at the tail between the port and the head
tailAngleColLog = strcmp(orientMatrixColIDs, 'TailAngle');
tailAngleBins = 0:15:180;
tailAngleOccupancyVector = histcounts(allOrientData(:,tailAngleColLog), tailAngleBins, 'Normalization', 'probability');
% Port Angle: Angle at the port between the head and the tail
portAngleColLog = strcmp(orientMatrixColIDs, 'PortAngle');
portAngleBins = 0:15:180;
portAngleOccupancyVector = histcounts(allOrientData(:,portAngleColLog), portAngleBins, 'Normalization', 'probability');

% Trial Position Occupancy Bins: Here I'm considering position to be a 2-D
% where the x-axis is Time and the Y-axis is trial position
posBins = 0.5:4.5;
timeBins = firingRateTrialData{1}(1,1)-timeBinSize:timeBinSize:firingRateTrialData{1}(end,1)+timeBinSize;
timeBinOccupancyMatrix = histcounts2(allTimeVals, allTrialIDvals, timeBins, posBins, 'Normalization', 'probability')';

%% Analyze the Orientation Data
% For this I will be breaking down the head and tail spatial position by
% trial position and analyzing the variation across and within trials. 
tPosHeadOccupyFig = PlotTrialPositionOccupancy(allOrientData, orientationTrialData, headXcolLog, headYcolLog, headXbins, headYbins, allTrialIDvals, odorIDsISC, timeBins, 'Head');
figure(tPosHeadOccupyFig);
annotation('textbox', [0.2 0.9 0.8 0.1], 'String', sprintf('Spatial Bin = %.02fmm, Temporal Bin = %.02fms', spatialBinSize*1.7, timeBinSize*1000), 'HorizontalAlignment', 'right', 'linestyle', 'none', 'FontSize', 12);
annotation('textbox', 'position', [0.01 0.01 0.9 0.05], 'string', sprintf('%s', cd), 'linestyle', 'none', 'interpreter', 'none');
orient(tPosHeadOccupyFig, 'tall');
orient(tPosHeadOccupyFig, 'landscape');
%     print
print('-painters', tPosHeadOccupyFig, '-dpdf', 'HeadLocationTrialPositionSummary');
tPosTailOccupyFig = PlotTrialPositionOccupancy(allOrientData, orientationTrialData, tailXcolLog, tailYcolLog, tailXbins, tailYbins, allTrialIDvals, odorIDsISC, timeBins, 'Tail');
figure(tPosTailOccupyFig);
annotation('textbox', [0.2 0.9 0.8 0.1], 'String', sprintf('Spatial Bin = %.02fmm, Temporal Bin = %.02fms', spatialBinSize*1.7, timeBinSize*1000), 'HorizontalAlignment', 'right', 'linestyle', 'none', 'FontSize', 12);
annotation('textbox', 'position', [0.01 0.01 0.9 0.05], 'string', sprintf('%s', cd), 'linestyle', 'none', 'interpreter', 'none');
orient(tPosTailOccupyFig, 'tall');
orient(tPosTailOccupyFig, 'landscape');
%     print
print('-painters', tPosTailOccupyFig, '-dpdf', 'TailLocationTrialPositionSummary');   
% Compute the deviation (distance) for the Head and Tail during the first
% frame of every trial from the spatial position on trial #1 in the
% sequence

% Compute the deviation (distance) of the Head and Tail during the trial
% period for every trial

% Calculate the spatial position for Head and Tail for each trial position
% for the entire trial period

% Calculate the difference across trial positions for the spatial position
% matrix (just do subtraction for now) and plot the difference in grayscale



%% Run Analyses
for uni = 1:length(ensembleUnitSummaries)
    %% Port Orientation Alone Analyses
    % Head Position Information
    [headFRmap, headICmap, headIC, headICrate] = CalculateFieldIC(allOrientData(:,headXcolLog), allOrientData(:,headYcolLog),...
        headXbins, headYbins, headOccupancyMatrix, allFRvals(:,uni), meanFRvals(uni));
    % Tail Position Information
    [tailFRmap, tailICmap, tailIC, tailICrate] = CalculateFieldIC(allOrientData(:,tailXcolLog), allOrientData(:,tailYcolLog),...
        tailXbins, tailYbins, tailOccupancyMatrix, allFRvals(:,uni), meanFRvals(uni));
    
    % Head Angle Information
    [headAngleFRvect, headAngleICvect, headAngleICval, headAngleICrate] = CalculateVectorIC(allOrientData(:,headAngleColLog),...
        headAngleBins, headAngleOccupancyVector, allFRvals(:,uni), meanFRvals(uni));
    % Tail Angle Information
    [tailAngleFRvect, tailAngleICvect, tailAngleICval, tailAngleICrate] = CalculateVectorIC(allOrientData(:,tailAngleColLog),...
        tailAngleBins, tailAngleOccupancyVector, allFRvals(:,uni), meanFRvals(uni));
    % Port Angle Information
    [portAngleFRvect, portAngleICvect, portAngleICval, portAngleICrate] = CalculateVectorIC(allOrientData(:,portAngleColLog),...
        portAngleBins, portAngleOccupancyVector, allFRvals(:,uni), meanFRvals(uni));
    
    
    figure
    % Head Occupancy Matrix
    subplot(3,3,1)
    imagesc(headOccupancyMatrix);
    set(gca, 'ydir', 'normal');
    title('Head Occupancy');
    cb = colorbar('southoutside');
    cb.Label.String = 'p(Occupancy)';
    z = colormap('jet');
    z(1,:) = [1 1 1];
    colormap(z);
    set(gca, 'clim', [-0.01, max(get(gca, 'clim'))]);
    % Head FR Map
    subplot(3,3,2)
    imagesc(headFRmap);
    set(gca, 'ydir', 'normal');
    title('Head Position FR');
    cb = colorbar('southoutside');
    cb.Label.String = 'Mean Firing Rate (spk/s)';
    colormap(z);
    % Head IC Map
    subplot(3,3,3)
    imagesc(headICmap);
    set(gca, 'ydir', 'normal');
    title(sprintf('Head Position IC (Overall = %.03f bits | %.03f bit/sec)', headIC, headICrate));
    cb = colorbar('southoutside');
    cb.Label.String = 'Information Content (bits)';
    colormap(z);
    
    % Tail Occupancy Matrix
    subplot(3,3,4)
    imagesc(tailOccupancyMatrix);
    set(gca, 'ydir', 'normal');
    title('Tail Occupancy');
    cb = colorbar('southoutside');
    cb.Label.String = 'p(Occupancy)';
    colormap(z);
    set(gca, 'clim', [-0.01, max(get(gca, 'clim'))]);
    % Head FR Map
    subplot(3,3,5)
    imagesc(tailFRmap);
    set(gca, 'ydir', 'normal');
    title('Tail Position FR');
    cb = colorbar('southoutside');
    cb.Label.String = 'Mean Firing Rate (spk/s)';
    colormap(z);
    % Head IC Map
    subplot(3,3,6)
    imagesc(tailICmap);
    set(gca, 'ydir', 'normal');
    title(sprintf('Tail Position IC (Overall = %.03f bits | %.03f bit/sec)', tailIC, tailICrate));
    cb = colorbar('southoutside');
    cb.Label.String = 'Information Content (bits)';
    colormap(z);
    
    % Head Angle
    subplot(3,3,7)
    yyaxis left
    plot(headAngleBins(2:end) - diff(headAngleBins)/2, headAngleOccupancyVector, 'linestyle', '--', 'color', 'k');
    hold on;
    yyaxis right
    plot(headAngleBins(2:end) - diff(headAngleBins)/2, headAngleFRvect, 'linestyle', '-', 'color', 'k');
    yyaxis left
    plot(headAngleBins(2:end) - diff(headAngleBins)/2, headAngleICvect, 'linestyle', '-', 'color', 'r');
    box off;
    yyaxis left
    ylabel 'Occupancy OR IC (bits)'
    yyaxis right
    ylabel 'Firing Rate'
    legend('Occ', 'FR', 'IC', 'location', 'best');
    title(sprintf('Head Angle (Overall IC = %.03f bits | %.03f bit/sec)', headAngleICval, headAngleICrate));
    
    subplot(3,3,8)
    yyaxis left
    plot(tailAngleBins(2:end) - diff(tailAngleBins)/2, tailAngleOccupancyVector, 'linestyle', '--', 'color', 'k');
    hold on;
    yyaxis right
    plot(tailAngleBins(2:end) - diff(tailAngleBins)/2, tailAngleFRvect, 'linestyle', '-', 'color', 'k');
    yyaxis left
    plot(tailAngleBins(2:end) - diff(tailAngleBins)/2, tailAngleICvect, 'linestyle', '-', 'color', 'r');
    box off;
    yyaxis left
    ylabel 'Occupancy OR IC (bits)'
    yyaxis right
    ylabel 'Firing Rate'
    legend('Occ', 'FR', 'IC', 'location', 'best');
    title(sprintf('Tail Angle (Overall IC = %.03f bits | %.03f bit/sec)', tailAngleICval, tailAngleICrate));

    subplot(3,3,9)
    yyaxis left
    plot(portAngleBins(2:end) - diff(portAngleBins)/2, portAngleOccupancyVector, 'linestyle', '--', 'color', 'k');
    hold on;
    yyaxis right
    plot(portAngleBins(2:end) - diff(portAngleBins)/2, portAngleFRvect, 'linestyle', '-', 'color', 'k');
    yyaxis left
    plot(portAngleBins(2:end) - diff(portAngleBins)/2, portAngleICvect, 'linestyle', '-', 'color', 'r');
    box off;
    yyaxis left
    ylabel 'Occupancy OR IC (bits)'
    yyaxis right
    ylabel 'Firing Rate'
    legend('Occ', 'FR', 'IC', 'location', 'best');
    title(sprintf('Port Angle (Overall IC = %.03f bits | %.03f bit/sec)', portAngleICval, portAngleICrate));
    
    annotation('textbox', [0.05 0.9 0.9 0.1], 'String', ensembleMatrixColIDs{uni+1}, 'linestyle', 'none', 'FontSize', 20); %#ok<USENS>
    annotation('textbox', [0.2 0.9 0.8 0.1], 'String', sprintf('Spatial Bin = %.02fmm, Gaussian = %.02fms', spatialBinSize*1.7, slideWindowSize), 'HorizontalAlignment', 'right', 'linestyle', 'none', 'FontSize', 12);
    annotation('textbox', 'position', [0.01 0.01 0.9 0.05], 'string', sprintf('%s', cd), 'linestyle', 'none', 'interpreter', 'none');
    drawnow
    
    orient(gcf, 'tall');
    orient(gcf, 'landscape');
%     print
    print('-painters', gcf, '-dpdf', sprintf('%s_OrientInfo_Summary', ensembleMatrixColIDs{uni+1}));
    close(gcf);
    
    %% Trial Position Alone
    % Trial Position
    [trialFRmap, trialICmap, trialIC, trialICrate] = CalculateFieldIC(allTimeVals, allTrialIDvals,...
        timeBins, posBins, timeBinOccupancyMatrix, allFRvals(:,uni), meanFRvals(uni));
    
    figure;
    subplot(1,3,1)
    imagesc(firingRateTrialData{1}(1,1):0.1:firingRateTrialData{1}(end,1),1:4,timeBinOccupancyMatrix);
    set(gca, 'ytick', 1:4);
    hold on;
    line([0 0], [0.5 4.5], 'linewidth', 2, 'color', 'k');
    title('Trial Occupancy');
    cb = colorbar('southoutside');
    cb.Label.String = 'p(Occupancy)';
    z = colormap('jet');
    z(1,:) = [1 1 1];
    colormap(z);
    set(gca, 'clim', [-0.01, max(get(gca, 'clim'))]);
    
    subplot(1,3,2)
    imagesc(firingRateTrialData{1}(1,1):0.1:firingRateTrialData{1}(end,1),1:4,trialFRmap);
    set(gca, 'ytick', 1:4);
    hold on;
    line([0 0], [0.5 4.5], 'linewidth', 2, 'color', 'k');
%     set(gca, 'ydir', 'normal');
    title('Trial Position FR');
    cb = colorbar('southoutside');
    cb.Label.String = 'Mean Firing Rate (spk/s)';
    colormap(z);
     
    subplot(1,3,3)
    imagesc(firingRateTrialData{1}(1,1):0.1:firingRateTrialData{1}(end,1),1:4,trialICmap);
    set(gca, 'ytick', 1:4);
    hold on;
    line([0 0], [0.5 4.5], 'linewidth', 2, 'color', 'k');
%     set(gca, 'ydir', 'normal');
    title(sprintf('Trial Position IC (Overall = %.03f bits | %.03f bit/sec)', trialIC, trialICrate));
    cb = colorbar('southoutside');
    cb.Label.String = 'Information Content (bits)';
    colormap(z);
    
    annotation('textbox', [0.05 0.9 0.9 0.1], 'String', ensembleMatrixColIDs{uni+1}, 'linestyle', 'none', 'FontSize', 20);
    annotation('textbox', [0.2 0.9 0.8 0.1], 'String', sprintf('Temporal Bin = %.02fmm, Gaussian = %.02fms', timeBinSize*1000, slideWindowSize), 'HorizontalAlignment', 'right', 'linestyle', 'none', 'FontSize', 12);
    annotation('textbox', 'position', [0.01 0.01 0.9 0.05], 'string', sprintf('%s', cd), 'linestyle', 'none', 'interpreter', 'none');
    drawnow
    
    orient(gcf, 'tall');
    orient(gcf, 'landscape');
%     print
    print('-painters', gcf, '-dpdf', sprintf('%s_PositionInfo_Summary', ensembleMatrixColIDs{uni+1}));
    close(gcf);
end

end
%%
function [frMap, icMap, totalICval, overallICrate] = CalculateFieldIC(curXdata, curYdata, xBinLims, yBinLims, occupancyMatrix, uniFR, meanFR)
frMap = nan(length(yBinLims)-1, length(xBinLims)-1);
icMap = nan(length(yBinLims)-1, length(xBinLims)-1);
for r = 1:length(yBinLims)-1
    yValLog = (curYdata>=yBinLims(r) & curYdata<yBinLims(r+1));
    for c = 1:length(xBinLims)-1
        xValLog = (curXdata>=xBinLims(c) & curXdata<xBinLims(c+1));
        frMap(r,c) = mean(uniFR(yValLog & xValLog));
        icMap(r,c) = occupancyMatrix(r,c)*(frMap(r,c)/meanFR)*log2(frMap(r,c)/meanFR);
    end
end
totalICval = nansum(icMap(:));
overallICrate = totalICval*meanFR;
end
%%
function [frVect, icVect, totalICval, overallICrate] = CalculateVectorIC(curXdata, xBinLims, occupancyVector, uniFR, meanFR)
frVect = nan(1,length(xBinLims)-1);
icVect = nan(1,length(xBinLims)-1);
for c = 1:length(xBinLims)-1
    curBinLog = (curXdata>=xBinLims(c) & curXdata<xBinLims(c+1));
    frVect(c) = mean(uniFR(curBinLog));
    icVect(c) = occupancyVector(c)*(frVect(c)/meanFR)*log2(frVect(c)/meanFR);
end
totalICval = nansum(icVect);
overallICrate = totalICval*meanFR;
end

%%
function [curTrialOccupancyPlot] = PlotTrialPositionOccupancy(allOrientData, orientationTrialData, xColLog, yColLog, xBins, yBins, allTrialIDvals, odorIDsISC, timeBins, feature)

trialOccupancy = cell(1,4);
for pos = 1:4
    trialOccupancy{pos} = histcounts2(allOrientData(allTrialIDvals==pos,xColLog), allOrientData(allTrialIDvals==pos,yColLog), xBins, yBins, 'Normalization', 'probability')';
    trialOccupancy{pos}(trialOccupancy{pos}==0) = nan;
end

% This is going to be hard coded unfortunately, but I think that's going to
% be the easiest way to handle it.
curTrialOccupancyPlot = figure;
jMod = colormap('jet');
jMod(1,:) = [1 1 1];
pMod = colormap('parula');
pMod(1,:) = [1 1 1];

pos1 = subplot(4,4,1);
imagesc(trialOccupancy{1});
set(pos1, 'ydir', 'normal', 'clim', [-0.01, max(get(gca, 'clim'))]);
title(sprintf('Trial 1 %s Location', feature));
colormap(pos1, jMod);

pos2 = subplot(4,4,6);
imagesc(trialOccupancy{2});
set(pos2, 'ydir', 'normal', 'clim', [-0.01, max(get(gca, 'clim'))]);
title(sprintf('Trial 2 %s Location', feature));
colormap(pos2, jMod);

pos3 = subplot(4,4,11);
imagesc(trialOccupancy{3});
set(pos3, 'ydir', 'normal', 'clim', [-0.01, max(get(gca, 'clim'))]);
title(sprintf('Trial 3 %s Location', feature));
colormap(pos3, jMod);

pos4 = subplot(4,4,16);
imagesc(trialOccupancy{4});
set(pos4, 'ydir', 'normal', 'clim', [-0.01, max(get(gca, 'clim'))]);
title(sprintf('Trial 4 %s Location', feature));
colormap(pos4, jMod);

% Now, zero out the nans... I know I nan'd out my zeros before. That's
% better for the bulk of the analyses, this is better for this analysis.
for t = 1:4
    trialOccupancy{t}(isnan(trialOccupancy{t})) = 0;
end

pos2m1 = subplot(4,4,2);
curDiffMtx = CalculateOccupyDiffMtx(trialOccupancy{1}, trialOccupancy{2});
imagesc(curDiffMtx);
set(pos2m1, 'ydir', 'normal', 'clim', [min(get(gca, 'clim'))-0.01, max(get(gca, 'clim'))]);
title('Trial 2 - Trial 1');
colormap(pos2m1, pMod);

pos3m1 = subplot(4,4,3);
curDiffMtx = CalculateOccupyDiffMtx(trialOccupancy{1}, trialOccupancy{3});
imagesc(curDiffMtx);
set(pos3m1, 'ydir', 'normal', 'clim', [min(get(gca, 'clim'))-0.01, max(get(gca, 'clim'))]);
title('Trial 3 - Trial 1');
colormap(pos3m1, pMod);

pos4m1 = subplot(4,4,4);
curDiffMtx = CalculateOccupyDiffMtx(trialOccupancy{1}, trialOccupancy{4});
imagesc(curDiffMtx);
set(pos4m1, 'ydir', 'normal', 'clim', [min(get(gca, 'clim'))-0.01, max(get(gca, 'clim'))]);
title('Trial 4 - Trial 1');
colormap(pos4m1, pMod);

pos3m2 = subplot(4,4,7);
curDiffMtx = CalculateOccupyDiffMtx(trialOccupancy{2}, trialOccupancy{3});
imagesc(curDiffMtx);
set(pos3m2, 'ydir', 'normal', 'clim', [min(get(gca, 'clim'))-0.01, max(get(gca, 'clim'))]);
title('Trial 3 - Trial 2');
colormap(pos3m2, pMod);

pos4m2 = subplot(4,4,8);
curDiffMtx = CalculateOccupyDiffMtx(trialOccupancy{2}, trialOccupancy{4});
imagesc(curDiffMtx);
set(pos4m2, 'ydir', 'normal', 'clim', [min(get(gca, 'clim'))-0.01, max(get(gca, 'clim'))]);
title('Trial 4 - Trial 2');
colormap(pos4m2, pMod);

pos4m3 = subplot(4,4,12);
curDiffMtx = CalculateOccupyDiffMtx(trialOccupancy{3}, trialOccupancy{4});
imagesc(curDiffMtx);
set(pos4m3, 'ydir', 'normal', 'clim', [min(get(gca, 'clim'))-0.01, max(get(gca, 'clim'))]);
title('Trial 4 - Trial 3');
colormap(pos4m3, pMod);

%% Extract spatial position in the first frame after poke in for every trial
xFF = nan(1,length(orientationTrialData));
yFF = nan(1,length(orientationTrialData));
trialMoveVect = cell(1,length(orientationTrialData));
for trl = 1:length(orientationTrialData)
    curTrlData = orientationTrialData{trl};
    curTrlData(isnan(curTrlData(:,2)),:) = [];
    xFF(trl) = curTrlData(find(curTrlData(:,1)>=0, 1, 'first'), xColLog);
    yFF(trl) = curTrlData(find(curTrlData(:,1)>=0, 1, 'first'), yColLog);
    trialMoveVect{trl} = [curTrlData(:,1), cell2mat(arrayfun(@(a,b)sqrt(((a-xFF(trl))^2)+((b-yFF(trl))^2)), curTrlData(:,xColLog), curTrlData(:,yColLog), 'uniformoutput', 0))];
end
% Compute the average spatial position for Head and Tail across trial
% positions
trialHeadXff = nan(2,4);
trialHeadYff = nan(2,4);
trialMoveVals = repmat({nan(2,length(timeBins)-1)}, [1 4]);
for pos = 1:4
    trialHeadXff(1,pos) = mean(xFF(odorIDsISC==pos));
    trialHeadXff(2,pos) = std(xFF(odorIDsISC==pos));
    
    trialHeadYff(1,pos) = mean(yFF(odorIDsISC==pos));
    trialHeadYff(2,pos) = std(yFF(odorIDsISC==pos));
    curPosTrlVects = cell2mat(trialMoveVect(odorIDsISC==pos)');
    tempDistVect = trialMoveVals{pos};
    for tb = 1:length(timeBins)-1
        tempDistVect(1,tb) = mean(curPosTrlVects(curPosTrlVects(:,1)>=timeBins(tb) & curPosTrlVects(:,1)<timeBins(tb+1),2));
        tempDistVect(2,tb) = std(curPosTrlVects(curPosTrlVects(:,1)>=timeBins(tb) & curPosTrlVects(:,1)<timeBins(tb+1),2));        
    end
    trialMoveVals{pos} = tempDistVect;
end

subplot(4,4,9)
errorbar(trialHeadXff(1,:), trialHeadYff(1,:),trialHeadYff(2,:),trialHeadYff(2,:), trialHeadXff(2,:), trialHeadXff(2,:), 'CapSize', 0);
for pos = 1:4
    text(trialHeadXff(1,pos), trialHeadYff(1,pos), sprintf('%i', pos), 'horizontalalignment', 'center', 'FontSize', 15, 'FontWeight', 'Bold');
end
ylabel('Y-Coordinate');
xlabel('X-Coordinate');
title(sprintf('Initial %s Position Upon Poke In', feature));

subplot(4,4,[13 14])
errorbar(timeBins(2:end)-mode(diff(timeBins)),trialMoveVals{1}(1,:), trialMoveVals{1}(2,:), 'capsize', 0, 'color', [44/255 168/255 224/255]);
hold on;
errorbar(timeBins(2:end)-mode(diff(timeBins)),trialMoveVals{2}(1,:), trialMoveVals{2}(2,:), 'capsize', 0, 'color', [154/255 133/255 122/255]);
errorbar(timeBins(2:end)-mode(diff(timeBins)),trialMoveVals{3}(1,:), trialMoveVals{3}(2,:), 'capsize', 0, 'color', [9/255 161/255 74/255]);
errorbar(timeBins(2:end)-mode(diff(timeBins)),trialMoveVals{4}(1,:), trialMoveVals{4}(2,:), 'capsize', 0, 'color', [128/255 66/255 151/255]);
legend('Position 1', 'Position 2', 'Position 3', 'Position 4', 'location', 'best');
set(gca, 'ylim', [-1 max(get(gca, 'ylim'))]);
title('Distance From First Poke In Frame');
xlabel 'Time'
ylabel 'Distance (pixels)'

    function diffOccupyMatrix = CalculateOccupyDiffMtx(occupyMatrix1, occupyMatrix2)
        diffOccupyMatrix = occupyMatrix2-occupyMatrix1;
        noOccupyMask = (occupyMatrix2==0 & occupyMatrix1==0);
        diffOccupyMatrix(noOccupyMask) = nan;
    end
end