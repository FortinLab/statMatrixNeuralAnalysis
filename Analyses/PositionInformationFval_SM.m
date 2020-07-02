
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
slideWindowSize = 100;
spatialBinSize = 5;
numPerms = 100;

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

% Now redefine the column ID vectors
portXcolLog = strcmp(orientMatrixColIDs, 'PortX');
portYcolLog = strcmp(orientMatrixColIDs, 'PortY');
headXcolLog = strcmp(orientMatrixColIDs, 'HeadX');
headYcolLog = strcmp(orientMatrixColIDs, 'HeadY');
tailXcolLog = strcmp(orientMatrixColIDs, 'TailX');
tailYcolLog = strcmp(orientMatrixColIDs, 'TailY');
headAngleColLog = strcmp(orientMatrixColIDs, 'HeadAngle');
tailAngleColLog = strcmp(orientMatrixColIDs, 'TailAngle');
portAngleColLog = strcmp(orientMatrixColIDs, 'PortAngle');
htLengthColLog = strcmp(orientMatrixColIDs, 'HeadTailLength');
hpLengthColLog = strcmp(orientMatrixColIDs, 'HeadPortLength');
ptLengthColLog = strcmp(orientMatrixColIDs, 'PortTailLength');

% Determine spatial bins
headXbins = floor(min(orientMatrix(:,headXcolLog))/spatialBinSize)*spatialBinSize:spatialBinSize:ceil(max(orientMatrix(:,headXcolLog))/spatialBinSize)*spatialBinSize;
headYbins = floor(min(orientMatrix(:,headYcolLog))/spatialBinSize)*spatialBinSize:spatialBinSize:ceil(max(orientMatrix(:,headYcolLog))/spatialBinSize)*spatialBinSize;

tailXbins = floor(min(orientMatrix(:,tailXcolLog))/spatialBinSize)*spatialBinSize:spatialBinSize:ceil(max(orientMatrix(:,tailXcolLog))/spatialBinSize)*spatialBinSize;
tailYbins = floor(min(orientMatrix(:,tailYcolLog))/spatialBinSize)*spatialBinSize:spatialBinSize:ceil(max(orientMatrix(:,tailYcolLog))/spatialBinSize)*spatialBinSize;

% Now Interpolate the port/head/tail locations. No interpolation of the
% other vals right now...
orientMatrix(:,portPosX) = interp1(orientMatrix(posIndices,1), orientMatrix(posIndices,portPosX), orientMatrix(:,1));
orientMatrix(:,portPosY) = interp1(orientMatrix(posIndices,1), orientMatrix(posIndices,portPosY), orientMatrix(:,1));
orientMatrix(:,headPosX) = interp1(orientMatrix(posIndices,1), orientMatrix(posIndices,headPosX), orientMatrix(:,1));
orientMatrix(:,headPosY) = interp1(orientMatrix(posIndices,1), orientMatrix(posIndices,headPosY), orientMatrix(:,1));
orientMatrix(:,tailPosX) = interp1(orientMatrix(posIndices,1), orientMatrix(posIndices,tailPosX), orientMatrix(:,1));
orientMatrix(:,tailPosY) = interp1(orientMatrix(posIndices,1), orientMatrix(posIndices,tailPosY), orientMatrix(:,1));

%% Calculate Instantaneous Firing Rates
% For this I'm opting to use the instantaneous firing rate (i.e. convolved
% with a gaussian) because the camera capture rate is not exactly uniform,
% which is what is needed for using discrete time bins. This may lead to
% some distortion in the spiking value... I may want to do it with discrete
% bins down the road to verify I get similar values.

% Convolve w/Gaussian
% Create Gaussian
instFRgauss = gausswin(slideWindowSize);
instFRgauss = instFRgauss/(length(instFRgauss)*mode(diff(behavMatrix(:,1))));

% Convolve Gaussian with Binary Spike Trains
uniInstFR = nan(size(ensembleMatrix,1), size(ensembleMatrix,2)-1);
for uni = 2:size(ensembleMatrix,2)
    uniInstFR(:,uni-1) = conv(ensembleMatrix(:,uni), instFRgauss, 'same');
end
% 
% % Alternatively, just bin the firing rates
% uniInstFR = [ensembleMatrix(:,1), nan(size(ensembleMatrix,1), size(ensembleMatrix,2)-1)];
% for uni = 2:size(ensembleMatrix,2)
%     uniInstFR(:,uni) = conv(ensembleMatrix(:,uni), ones(1,slideWindowSize)/(length(instFRgauss)*mode(diff(behavMatrix(:,1)))), 'same'); 
% end

%% Organize Trial Data
% Create Extraction Matrices
pokeInTrialMatrix = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [-0.5 0.5], 'PokeIn');
pokeInTSs = behavMatrix(pokeInTrialMatrix(1).TrialLogVect,1)-behavMatrix(pokeInTrialMatrix(1).PokeInIndex,1);
pokeOutTrialMatrix = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [-0.5 0.5], 'PokeOut');
pokeOutTSs = behavMatrix(pokeOutTrialMatrix(1).TrialLogVect,1)-behavMatrix(pokeOutTrialMatrix(1).PokeOutIndex,1);

% Create Trial Logical Vectors
corrTrlLog = [pokeInTrialMatrix.Performance];
isLog = [pokeInTrialMatrix.TranspositionDistance]==0;
pos1log = [pokeInTrialMatrix.Position]==1;
trialLog = (corrTrlLog & isLog & ~pos1log);
seqType = 'All InSeq Correct Except Position 1';
% trialLog = (corrTrlLog & isLog);
% seqType = 'All InSeq Correct Trials';
% trialLog = (corrTrlLog);
% seqType = 'All Correct Trials';
posIDsISC = [pokeInTrialMatrix(trialLog).Position];

% Pull out spiking activity
pokeInFR = ExtractTrialData_SM(pokeInTrialMatrix, uniInstFR); %#ok<*NODEF>
pokeInFRisc = cell2mat(reshape(pokeInFR(trialLog), [1,1,sum(trialLog)]));
pokeOutFR = ExtractTrialData_SM(pokeOutTrialMatrix, uniInstFR);
pokeOutFRisc = cell2mat(reshape(pokeOutFR(trialLog), [1,1,sum(trialLog)]));

% Pull out and interpolate positional information
pokeInOrientation = ExtractTrialData_SM(pokeInTrialMatrix, orientMatrix); %#ok<*NODEF>
pokeInOrientationInterp = cell(size(pokeInOrientation));
pokeOutOrientation = ExtractTrialData_SM(pokeOutTrialMatrix, orientMatrix);
pokeOutOrientationInterp = cell(size(pokeOutOrientation));

%% Remove Units with low firing rates
uniFRthreshLog = max([max(mean(pokeInFRisc,3)); max(mean(pokeOutFRisc,3))])<1;
pokeInFRisc(:,uniFRthreshLog,:) = [];
pokeOutFRisc(:,uniFRthreshLog,:) = [];
goodUniNames = {ensembleUnitSummaries(~uniFRthreshLog).UnitName};

%% Pre-Process the positional data v2
% First thing to do is nan out the pre trial periods that don't have
% positional data (in the pokeIn aligned) as well as the post trial periods
% that don't have positional data (in the pokeOut aligned)

for trl = 1:length(pokeInOrientation)
    curPIori = pokeInOrientation{trl};
    firstPreTrialNdx = find(~isnan(curPIori(:,headAngleColLog)),1,'first');
    curPIori(1:firstPreTrialNdx, 2:7) = nan;
    for ndx = firstPreTrialNdx:size(curPIori,1)
        curPIori(ndx,ptLengthColLog) = sqrt((curPIori(ndx,tailXcolLog) - curPIori(ndx,portXcolLog))^2 + (curPIori(ndx,tailYcolLog) - curPIori(ndx,portYcolLog))^2);
        curPIori(ndx,hpLengthColLog) = sqrt((curPIori(ndx,portXcolLog) - curPIori(ndx,headXcolLog))^2 + (curPIori(ndx,portYcolLog) - curPIori(ndx,headYcolLog))^2);
        curPIori(ndx,htLengthColLog) = sqrt((curPIori(ndx,tailXcolLog) - curPIori(ndx,headXcolLog))^2 + (curPIori(ndx,tailYcolLog) - curPIori(ndx,headYcolLog))^2);
        curPIori(ndx,headAngleColLog) = rad2deg(acos((curPIori(ndx,htLengthColLog)^2 + curPIori(ndx,hpLengthColLog)^2 - curPIori(ndx,ptLengthColLog)^2)/(2*curPIori(ndx,htLengthColLog)*curPIori(ndx,hpLengthColLog))));
        curPIori(ndx,portAngleColLog) = rad2deg(acos((curPIori(ndx,ptLengthColLog)^2 + curPIori(ndx,hpLengthColLog)^2 - curPIori(ndx,htLengthColLog)^2)/(2*curPIori(ndx,ptLengthColLog)*curPIori(ndx,hpLengthColLog))));
        curPIori(ndx,tailAngleColLog) = rad2deg(acos((curPIori(ndx,htLengthColLog)^2 + curPIori(ndx,ptLengthColLog)^2 - curPIori(ndx,hpLengthColLog)^2)/(2*curPIori(ndx,htLengthColLog)*curPIori(ndx,ptLengthColLog)))); 
    end
    curPOori = pokeOutOrientation{trl};
    lastTrialNdx = find(~isnan(curPOori(:,headAngleColLog)),1,'last');
    curPOori(lastTrialNdx+1:end, 2:7) = nan;
    for ndx = 1:lastTrialNdx
        curPOori(ndx,ptLengthColLog) = sqrt((curPOori(ndx,tailXcolLog) - curPOori(ndx,portXcolLog))^2 + (curPOori(ndx,tailYcolLog) - curPOori(ndx,portYcolLog))^2);
        curPOori(ndx,hpLengthColLog) = sqrt((curPOori(ndx,portXcolLog) - curPOori(ndx,headXcolLog))^2 + (curPOori(ndx,portYcolLog) - curPOori(ndx,headYcolLog))^2);
        curPOori(ndx,htLengthColLog) = sqrt((curPOori(ndx,tailXcolLog) - curPOori(ndx,headXcolLog))^2 + (curPOori(ndx,tailYcolLog) - curPOori(ndx,headYcolLog))^2);
        curPOori(ndx,headAngleColLog) = rad2deg(acos((curPOori(ndx,htLengthColLog)^2 + curPOori(ndx,hpLengthColLog)^2 - curPOori(ndx,ptLengthColLog)^2)/(2*curPOori(ndx,htLengthColLog)*curPOori(ndx,hpLengthColLog))));
        curPOori(ndx,portAngleColLog) = rad2deg(acos((curPOori(ndx,ptLengthColLog)^2 + curPOori(ndx,hpLengthColLog)^2 - curPOori(ndx,htLengthColLog)^2)/(2*curPOori(ndx,ptLengthColLog)*curPOori(ndx,hpLengthColLog))));
        curPOori(ndx,tailAngleColLog) = rad2deg(acos((curPOori(ndx,htLengthColLog)^2 + curPOori(ndx,ptLengthColLog)^2 - curPOori(ndx,hpLengthColLog)^2)/(2*curPOori(ndx,htLengthColLog)*curPOori(ndx,ptLengthColLog)))); 
    end
    
    htVal(curPosNdx) = sqrt((curTailX - curHeadX)^2 + (curTailY - curHeadY)^2);
    hpVal(curPosNdx) = sqrt((curPortX - curHeadX)^2 + (curPortY - curHeadY)^2);
    ptVal(curPosNdx) = sqrt((curTailX - curPortX)^2 + (curTailY - curPortY)^2);

    pokeInOrientationInterp{trl} = curPIori;
    pokeOutOrientationInterp{trl} = curPOori;
end
% Now bin the orientation data
orientBinValsPI = cell(size(pokeInOrientationInterp));
orientBinValsPO = cell(size(pokeOutOrientationInterp));
for trl = 1:length(pokeInOrientation)
    curPIori = pokeInOrientationInterp{trl};
    tempPIndxVals = nan(size(curPIori,1),2);
    for piT = 1:size(curPIori,1)
        if ~isnan(curPIori(piT,headXcolLog))
            tempPIndxVals(piT,1) = sub2ind([length(headYbins)-1 length(headXbins)-1], find(curPIori(piT,headYcolLog)>=headYbins,1,'last'), find(curPIori(piT,headXcolLog)>=headXbins,1,'last'));
            tempPIndxVals(piT,2) = sub2ind([length(tailYbins)-1 length(tailXbins)-1], find(curPIori(piT,tailYcolLog)>=tailYbins,1,'last'), find(curPIori(piT,tailXcolLog)>=tailXbins,1,'last'));
        end
    end
    curPOori = pokeOutOrientationInterp{trl};
    tempPOndxVals = nan(size(curPOori,1),2);
    for poT = 1:size(curPOori,1)
        if ~isnan(curPOori(poT,headXcolLog))
            tempPOndxVals(poT,1) = sub2ind([length(headYbins)-1 length(headXbins)-1], find(curPOori(poT,headYcolLog)>=headYbins,1,'last'), find(curPOori(poT,headXcolLog)>=headXbins,1,'last'));
            tempPOndxVals(poT,2) = sub2ind([length(tailYbins)-1 length(tailXbins)-1], find(curPOori(poT,tailYcolLog)>=tailYbins,1,'last'), find(curPOori(poT,tailXcolLog)>=tailXbins,1,'last'));
        end
    end
    orientBinValsPI{trl} = tempPIndxVals;
    orientBinValsPO{trl} = tempPOndxVals;
end    
pokeInOrientISC = cell2mat(reshape(pokeInOrientationInterp(trialLog), [1,1,sum(trialLog)]));
pokeOutOrientISC = cell2mat(reshape(pokeOutOrientationInterp(trialLog), [1,1,sum(trialLog)]));

pokeInOrientBinsISC = cell2mat(reshape(orientBinValsPI(trialLog), [1,1,sum(trialLog)]));
pokeOutOrientBinsISC = cell2mat(reshape(orientBinValsPO(trialLog), [1,1,sum(trialLog)]));

trialOrientBinsISC = [pokeInOrientBinsISC; pokeOutOrientBinsISC];

fValsPI = cell(4,length(ensembleUnitSummaries));
fValsPO = cell(4,length(ensembleUnitSummaries));

%% Calculate each unit's FR map and positional firing rates and such...
for uni = 1:length(goodUniNames)
    curUniFR = [pokeInFRisc(:,uni,:); pokeOutFRisc(:,uni,:)];
    tempHeadFRmap = nan(length(headYbins)-1, length(headXbins)-1);
    for headNdx = 1:((length(headXbins)-1)*(length(headYbins)-1))
        tempHeadFRmap(headNdx) = mean(curUniFR(trialOrientBinsISC(:,1,:)==headNdx));
    end
    tempTailFRmap = nan(length(tailYbins)-1, length(tailXbins)-1);
    for tailNdx = 1:((length(tailXbins)-1)*(length(tailYbins)-1))
        tempTailFRmap(tailNdx) = mean(curUniFR(trialOrientBinsISC(:,2,:)==tailNdx));
    end
    
    figure;
    sp1 = subplot(1,2,1);
    imagesc(tempHeadFRmap)
    set(gca, 'ydir', 'normal');
    title('Head Firing Rate Map');
    cb = colorbar('southoutside');
    cb.Label.String = 'Mean Firing Rate (spk/s)';
    z = colormap('jet');
    z(1,:) = [1 1 1];
    colormap(z);
    set(gca, 'clim', [-0.01, max(get(gca, 'clim'))]);
    sp2 = subplot(1,2,2);
    imagesc(tempTailFRmap);
    set(gca, 'ydir', 'normal');
    title('Tail Firing Rate Map');
    cb = colorbar('southoutside');
    cb.Label.String = 'Mean Firing Rate (spk/s)';
    colormap(z);
    set(sp1, 'clim', [-0.01, max([get(sp1, 'clim'), get(sp2, 'clim')])]);
    set(sp2, 'clim', [-0.01, max([get(sp1, 'clim'), get(sp2, 'clim')])]);
    annotation('textbox', [0.05 0.9 0.9 0.1], 'String', goodUniNames{uni}, 'linestyle', 'none', 'FontSize', 20);
    annotation('textbox', [0.2 0.9 0.8 0.1], 'String', sprintf('%s Spatial Bin = %.02fmm Gaussian = %.02fms Number of Permutations = %i', seqType, spatialBinSize*1.7, slideWindowSize, numPerms), 'HorizontalAlignment', 'right', 'linestyle', 'none', 'FontSize', 12);
    annotation('textbox', 'position', [0.01 0.01 0.9 0.05], 'string', sprintf('%s', cd), 'linestyle', 'none', 'interpreter', 'none');
    
    orient(gcf, 'tall');
    orient(gcf, 'landscape');
    print('-painters', gcf, '-dpdf', sprintf('%s_FiringRateLocationMaps', goodUniNames{uni}));

    %*****************Make these figures broken down by sequence position
    % Then do a comparison of FR maps across sequence position to look for
    % modifications in FR maps by sequence position.
end

%% Run Sliding F-Ratio Analysis
for uni = 1:length(goodUniNames)
    curUniPokeInFR = reshape(pokeInFRisc(:,uni,:), [size(pokeInFRisc,1), size(pokeInFRisc,3), 1]);
    curUniPokeOutFR = reshape(pokeOutFRisc(:,uni,:), [size(pokeOutFRisc,1), size(pokeOutFRisc,3), 1]);
    
    % Poke In Oriented
    [curUniPIfvalPOS, curUniPIfvalPOSz] = SlidingFvalCalc(curUniPokeInFR, posIDsISC, numPerms);
    fValsPI{1,uni} = curUniPIfvalPOSz;
    [curUniPIfvalHeadLoc, curUniPIfvalHeadLocZ] = SlidingFvalCalc(curUniPokeInFR, reshape(pokeInOrientBinsISC(:,1,:), [size(pokeInOrientBinsISC,1), size(pokeInOrientBinsISC,3), 1]), numPerms);
    fValsPI{2,uni} = curUniPIfvalHeadLocZ;
    [curUniPIfvalTailLoc, curUniPIfvalTailLocZ] = SlidingFvalCalc(curUniPokeInFR, reshape(pokeInOrientBinsISC(:,2,:), [size(pokeInOrientBinsISC,1), size(pokeInOrientBinsISC,3), 1]), numPerms);
    fValsPI{3,uni} = curUniPIfvalTailLocZ;
    
    [curUniPIfvalInteract, curUniPIfvalInteractZ] = SlidingFvalCalcInteract(curUniPokeInFR, posIDsISC,...
        reshape(pokeInOrientBinsISC(:,1,:), [size(pokeInOrientBinsISC,1), size(pokeInOrientBinsISC,3), 1]),...
        reshape(pokeInOrientBinsISC(:,2,:), [size(pokeInOrientBinsISC,1), size(pokeInOrientBinsISC,3), 1]),...
        numPerms);
    fValsPI{4,uni} = curUniPIfvalInteractZ;
    
    
    %
    %     [curUniPIfvalHeadX, curUniPIfvalHeadXz] = SlidingFvalCalc(curUniPokeInFR, reshape(pokeInOrientISC(:,headXcolLog,:), [size(pokeInOrientISC,1), size(pokeInOrientISC,3), 1]), numPerms);
    %     [curUniPIfvalHeadY, curUniPIfvalHeadYz] = SlidingFvalCalc(curUniPokeInFR, reshape(pokeInOrientISC(:,headYcolLog,:), [size(pokeInOrientISC,1), size(pokeInOrientISC,3), 1]), numPerms);
    %
    %     [curUniPIfvalTailX, curUniPIfvalTailXz] = SlidingFvalCalc(curUniPokeInFR, reshape(pokeInOrientISC(:,tailXcolLog,:), [size(pokeInOrientISC,1), size(pokeInOrientISC,3), 1]), numPerms);
    %     [curUniPIfvalTailY, curUniPIfvalTailYz] = SlidingFvalCalc(curUniPokeInFR, reshape(pokeInOrientISC(:,tailYcolLog,:), [size(pokeInOrientISC,1), size(pokeInOrientISC,3), 1]), numPerms);
    %
    %     [curUniPIfvalPortAngle, curUniPIfvalPortAnglez] = SlidingFvalCalc(curUniPokeInFR, reshape(pokeInOrientISC(:,portAngleColLog,:), [size(pokeInOrientISC,1), size(pokeInOrientISC,3), 1]), numPerms);
    %     [curUniPIfvalHeadAngle, curUniPIfvalHeadAnglez] = SlidingFvalCalc(curUniPokeInFR, reshape(pokeInOrientISC(:,headAngleColLog,:), [size(pokeInOrientISC,1), size(pokeInOrientISC,3), 1]), numPerms);
    %     [curUniPIfvalTailAngle, curUniPIfvalTailAnglez] = SlidingFvalCalc(curUniPokeInFR, reshape(pokeInOrientISC(:,tailAngleColLog,:), [size(pokeInOrientISC,1), size(pokeInOrientISC,3), 1]), numPerms);
    %
    %     [curUniPIfvalHeadTailDist, curUniPIfvalHeadTailDistZ] = SlidingFvalCalc(curUniPokeInFR, reshape(pokeInOrientISC(:,htLengthColLog,:), [size(pokeInOrientISC,1), size(pokeInOrientISC,3), 1]), numPerms);
    %     [curUniPIfvalHeadPortDist, curUniPIfvalHeadPortDistZ] = SlidingFvalCalc(curUniPokeInFR, reshape(pokeInOrientISC(:,hpLengthColLog,:), [size(pokeInOrientISC,1), size(pokeInOrientISC,3), 1]), numPerms);
    %     [curUniPIfvalTailPortDist, curUniPIfvalTailPortDistZ] = SlidingFvalCalc(curUniPokeInFR, reshape(pokeInOrientISC(:,ptLengthColLog,:), [size(pokeInOrientISC,1), size(pokeInOrientISC,3), 1]), numPerms);
    %
    % Poke Out Oriented
    [curUniPOfvalPOS, curUniPOfvalPOSz] = SlidingFvalCalc(curUniPokeOutFR, posIDsISC, numPerms);
    fValsPO{1,uni} = curUniPOfvalPOSz;
    [curUniPOfvalHeadLoc, curUniPOfvalHeadLocZ] = SlidingFvalCalc(curUniPokeOutFR, reshape(pokeOutOrientBinsISC(:,1,:), [size(pokeOutOrientBinsISC,1), size(pokeOutOrientBinsISC,3), 1]), numPerms);
    fValsPO{2,uni} = curUniPOfvalHeadLocZ;
    [curUniPOfvalTailLoc, curUniPOfvalTailLocZ] = SlidingFvalCalc(curUniPokeOutFR, reshape(pokeOutOrientBinsISC(:,2,:), [size(pokeOutOrientBinsISC,1), size(pokeOutOrientBinsISC,3), 1]), numPerms);
    fValsPO{3,uni} = curUniPOfvalTailLocZ;
    
    [curUniPOfvalInteract, curUniPOfvalInteractZ] = SlidingFvalCalcInteract(curUniPokeOutFR, posIDsISC,...
        reshape(pokeOutOrientBinsISC(:,1,:), [size(pokeOutOrientBinsISC,1), size(pokeOutOrientBinsISC,3), 1]),...
        reshape(pokeOutOrientBinsISC(:,2,:), [size(pokeOutOrientBinsISC,1), size(pokeOutOrientBinsISC,3), 1]),...
        numPerms);
    fValsPO{4,uni} = curUniPOfvalInteractZ;
    %
    %     [curUniPOfvalHeadX, curUniPOfvalHeadXz] = SlidingFvalCalc(curUniPokeOutFR, reshape(pokeOutOrientISC(:,headXcolLog,:), [size(pokeOutOrientISC,1), size(pokeOutOrientISC,3), 1]), numPerms);
    %     [curUniPOfvalHeadY, curUniPOfvalHeadYz] = SlidingFvalCalc(curUniPokeOutFR, reshape(pokeOutOrientISC(:,headYcolLog,:), [size(pokeOutOrientISC,1), size(pokeOutOrientISC,3), 1]), numPerms);
    %
    %     [curUniPOfvalTailX, curUniPOfvalTailXz] = SlidingFvalCalc(curUniPokeOutFR, reshape(pokeOutOrientISC(:,tailXcolLog,:), [size(pokeOutOrientISC,1), size(pokeOutOrientISC,3), 1]), numPerms);
    %     [curUniPOfvalTailY, curUniPOfvalTailYz] = SlidingFvalCalc(curUniPokeOutFR, reshape(pokeOutOrientISC(:,tailYcolLog,:), [size(pokeOutOrientISC,1), size(pokeOutOrientISC,3), 1]), numPerms);
    %
    %     [curUniPOvalPortAngle, curUniPOfvalPortAnglez] = SlidingFvalCalc(curUniPokeOutFR, reshape(pokeOutOrientISC(:,portAngleColLog,:), [size(pokeOutOrientISC,1), size(pokeOutOrientISC,3), 1]), numPerms);
    %     [curUniPOfvalHeadAngle, curUniPOfvalHeadAnglez] = SlidingFvalCalc(curUniPokeOutFR, reshape(pokeOutOrientISC(:,headAngleColLog,:), [size(pokeOutOrientISC,1), size(pokeOutOrientISC,3), 1]), numPerms);
    %     [curUniPOfvalTailAngle, curUniPOfvalTailAnglez] = SlidingFvalCalc(curUniPokeOutFR, reshape(pokeOutOrientISC(:,tailAngleColLog,:), [size(pokeOutOrientISC,1), size(pokeOutOrientISC,3), 1]), numPerms);
    %
    %     [curUniPOfvalHeadTailDist, curUniPOfvalHeadTailDistZ] = SlidingFvalCalc(curUniPokeOutFR, reshape(pokeOutOrientISC(:,htLengthColLog,:), [size(pokeOutOrientISC,1), size(pokeOutOrientISC,3), 1]), numPerms);
    %     [curUniPOfvalHeadPortDist, curUniPOfvalHeadPortDistZ] = SlidingFvalCalc(curUniPokeOutFR, reshape(pokeOutOrientISC(:,hpLengthColLog,:), [size(pokeOutOrientISC,1), size(pokeOutOrientISC,3), 1]), numPerms);
    %     [curUniPOfvalTailPortDist, curUniPOfvalTailPortDistZ] = SlidingFvalCalc(curUniPokeOutFR, reshape(pokeOutOrientISC(:,ptLengthColLog,:), [size(pokeOutOrientISC,1), size(pokeOutOrientISC,3), 1]), numPerms);
    %
    figure;
    sp1 = subplot(2,1,1);
    plot(pokeInTSs, curUniPIfvalPOSz,'linewidth', 2);
    hold on;
    plot(pokeInTSs, curUniPIfvalHeadLocZ, 'k');
    plot(pokeInTSs, curUniPIfvalTailLocZ, '--k');
    plot(pokeInTSs, curUniPIfvalInteractZ, 'r');
    %     plot(pokeInTSs, curUniPIfvalHeadXz, 'k');
    %     plot(pokeInTSs, curUniPIfvalHeadYz, '--k');
    %     plot(pokeInTSs, curUniPIfvalTailXz, 'r');
    %     plot(pokeInTSs, curUniPIfvalTailYz, '--r');
    %     plot(pokeInTSs, curUniPIfvalPortAnglez, 'g');
    %     plot(pokeInTSs, curUniPIfvalHeadAnglez, '--g');
    %     plot(pokeInTSs, curUniPIfvalTailAnglez, '.g');
    %     plot(pokeInTSs, curUniPIfvalHeadTailDistZ, 'c');
    %     plot(pokeInTSs, curUniPIfvalHeadPortDistZ, '--c');
    %     plot(pokeInTSs, curUniPIfvalTailPortDistZ, '.c');
    title('Poke In Aligned')
    axis tight
    legend('Sequence Position', 'Head', 'Tail', 'Interaction', 'location', 'best');
    %     legend('Sequence Position', 'Head Loc X', 'Head Loc Y', 'Tail Loc X', 'Tail Loc Y', 'Port Angle', 'Head Angle', 'Tail Angle', 'Head-Tail', 'Head-Port', 'Tail-Port', 'location', 'best');
    
    sp2 = subplot(2,1,2);
    plot(pokeOutTSs, curUniPOfvalPOSz, 'linewidth', 2);
    hold on;
    plot(pokeOutTSs, curUniPOfvalHeadLocZ, 'k');
    plot(pokeOutTSs, curUniPOfvalTailLocZ, '--k');
    plot(pokeOutTSs, curUniPOfvalInteractZ, 'r');
    %     plot(pokeOutTSs, curUniPOfvalHeadXz, 'k');
    %     plot(pokeOutTSs, curUniPOfvalHeadYz, '--k');
    %     plot(pokeOutTSs, curUniPOfvalTailXz, 'r');
    %     plot(pokeOutTSs, curUniPOfvalTailYz, '--r');
    %     plot(pokeOutTSs, curUniPOfvalPortAnglez, 'g');
    %     plot(pokeOutTSs, curUniPOfvalHeadAnglez, '--g');
    %     plot(pokeOutTSs, curUniPOfvalTailAnglez, '.g');
    %     plot(pokeOutTSs, curUniPOfvalHeadTailDistZ, 'c');
    %     plot(pokeOutTSs, curUniPOfvalHeadPortDistZ, '--c');
    %     plot(pokeOutTSs, curUniPOfvalTailPortDistZ, '.c');
    title('Poke Out Aligned');
    axis tight
    
    linkaxes([sp1, sp2], 'y');
    annotation('textbox', [0.05 0.9 0.9 0.1], 'String', goodUniNames{uni}, 'linestyle', 'none', 'FontSize', 20);
    annotation('textbox', [0.2 0.9 0.8 0.1], 'String', sprintf('%s Spatial Bin = %.02fmm Gaussian = %.02fms Number of Permutations = %i', seqType, spatialBinSize*1.7, slideWindowSize, numPerms), 'HorizontalAlignment', 'right', 'linestyle', 'none', 'FontSize', 12);
    annotation('textbox', 'position', [0.01 0.01 0.9 0.05], 'string', sprintf('%s', cd), 'linestyle', 'none', 'interpreter', 'none');
    
    orient(gcf, 'tall');
    orient(gcf, 'landscape');
    %     print
    print('-painters', gcf, '-dpdf', sprintf('%s_OrientFvalInfo_Summary', goodUniNames{uni}));
end

%% Now summarize the population data
% Row 1 = Position vs Head Location
% Row 2 = Position vs Tail Location
% Row 3 = Head vs Tail Location
piCorr = nan(3,length(goodUniNames));
poCorr = nan(3,length(goodUniNames));
allCorr = nan(3,length(goodUniNames));
for uni = 1:length(goodUniNames)
    curPIcorr = corr(cell2mat(fValsPI(:,uni))');
    piCorr(1,uni) = curPIcorr(2,1);
    piCorr(2,uni) = curPIcorr(3,1);
    piCorr(3,uni) = curPIcorr(3,2);
    
    curPOcorr = corr(cell2mat(fValsPO(:,uni))');
    poCorr(1,uni) = curPOcorr(2,1);
    poCorr(2,uni) = curPOcorr(3,1);
    poCorr(3,uni) = curPOcorr(3,2);
    
    curAllCorr = corr([cell2mat(fValsPI(:,uni)), cell2mat(fValsPO(:,uni))]');
    allCorr(1,uni) = curAllCorr(2,1);
    allCorr(2,uni) = curAllCorr(3,1);
    allCorr(3,uni) = curAllCorr(3,2);
end
figure; 
sp1 = subplot(3,3,1);
histogram(piCorr(1,:), -1:0.1:1);
hold on;
line(gca, [nanmean(piCorr(1,:)) nanmean(piCorr(1,:))], get(gca, 'ylim'), 'color','r', 'linewidth', 2);
title('PI: Position vs. Head');
sp2 = subplot(3,3,2);
histogram(piCorr(2,:), -1:0.1:1);
hold on;
line(gca, [nanmean(piCorr(2,:)) nanmean(piCorr(2,:))], get(gca, 'ylim'), 'color','r', 'linewidth', 2);
title('PI: Position vs. Tail');
sp3 = subplot(3,3,3);
histogram(piCorr(3,:), -1:0.1:1);
hold on;
line(gca, [nanmean(piCorr(3,:)) nanmean(piCorr(3,:))], get(gca, 'ylim'), 'color','r', 'linewidth', 2);
title('PI: Head vs. Tail');
sp4 = subplot(3,3,4);
histogram(poCorr(1,:), -1:0.1:1);
hold on;
line(gca, [nanmean(poCorr(1,:)) nanmean(poCorr(1,:))], get(gca, 'ylim'), 'color','r', 'linewidth', 2);
title('PO: Position vs. Head');
sp5 = subplot(3,3,5);
histogram(poCorr(2,:), -1:0.1:1);
hold on;
line(gca, [nanmean(poCorr(2,:)) nanmean(poCorr(2,:))], get(gca, 'ylim'), 'color','r', 'linewidth', 2);
title('PO: Position vs. Tail');
sp6 = subplot(3,3,6);
histogram(poCorr(3,:), -1:0.1:1);
hold on;
line(gca, [nanmean(poCorr(3,:)) nanmean(poCorr(3,:))], get(gca, 'ylim'), 'color','r', 'linewidth', 2);
title('PO: Head vs. Tail');
sp7 = subplot(3,3,7);
histogram(allCorr(1,:), -1:0.1:1);
hold on;
line(gca, [nanmean(allCorr(1,:)) nanmean(allCorr(1,:))], get(gca, 'ylim'), 'color','r', 'linewidth', 2);
title('All: Position vs. Head');
sp8 = subplot(3,3,8);
histogram(allCorr(2,:), -1:0.1:1);
hold on;
line(gca, [nanmean(allCorr(2,:)) nanmean(allCorr(2,:))], get(gca, 'ylim'), 'color','r', 'linewidth', 2);
title('All: Position vs. Tail');
sp9 = subplot(3,3,9);
histogram(allCorr(3,:), -1:0.1:1);
hold on;
line(gca, [nanmean(allCorr(3,:)) nanmean(allCorr(3,:))], get(gca, 'ylim'), 'color','r', 'linewidth', 2);
title('All: Head vs. Tail');

linkaxes([sp1, sp2, sp2, sp3, sp4, sp5, sp6, sp7, sp8, sp9], 'y');

orient(gcf, 'tall');
orient(gcf, 'landscape');
annotation('textbox', [0.05 0.9 0.9 0.1], 'String', 'All Timepoints', 'linestyle', 'none', 'FontSize', 20);
annotation('textbox', [0.2 0.9 0.8 0.1], 'String', sprintf('%s Spatial Bin = %.02fmm Gaussian = %.02fms Number of Permutations = %i', seqType, spatialBinSize*1.7, slideWindowSize, numPerms), 'HorizontalAlignment', 'right', 'linestyle', 'none', 'FontSize', 12);
annotation('textbox', 'position', [0.01 0.01 0.9 0.05], 'string', sprintf('%s', cd), 'linestyle', 'none', 'interpreter', 'none');
print('-painters', gcf, '-dpdf', 'LocationVsPositionCorrelationSummary(All)');
%% Now run a thresholded version
% Row 1 = Position vs Head Location
% Row 2 = Position vs Tail Location
% Row 3 = Head vs Tail Location
piCorr = nan(3,length(goodUniNames));
poCorr = nan(3,length(goodUniNames));
allCorr = nan(3,length(goodUniNames));
for uni = 1:length(goodUniNames)
    curPIphCorr = corrcoef(fValsPI{1,uni}(fValsPI{1,uni}>1 | fValsPI{2,uni}>1), fValsPI{2,uni}(fValsPI{1,uni}>1 | fValsPI{2,uni}>1));
    piCorr(1,uni) = curPIphCorr(2,1);
    curPIptCorr = corrcoef(fValsPI{1,uni}(fValsPI{1,uni}>1 | fValsPI{3,uni}>1), fValsPI{3,uni}(fValsPI{1,uni}>1 | fValsPI{3,uni}>1));
    piCorr(2,uni) = curPIptCorr(2,1);
    curPIhtCorr = corrcoef(fValsPI{2,uni}(fValsPI{2,uni}>1 | fValsPI{3,uni}>1), fValsPI{3,uni}(fValsPI{2,uni}>1 | fValsPI{3,uni}>1));
    piCorr(3,uni) = curPIhtCorr(2,1);
    
    curPOphCorr = corrcoef(fValsPO{1,uni}(fValsPO{1,uni}>1 | fValsPO{2,uni}>1), fValsPO{2,uni}(fValsPO{1,uni}>1 | fValsPO{2,uni}>1));
    poCorr(1,uni) = curPOphCorr(2,1);
    curPOptCorr = corrcoef(fValsPO{1,uni}(fValsPO{1,uni}>1 | fValsPO{3,uni}>1), fValsPO{3,uni}(fValsPO{1,uni}>1 | fValsPO{3,uni}>1));
    poCorr(2,uni) = curPOptCorr(2,1);
    curPOhtCorr = corrcoef(fValsPO{2,uni}(fValsPO{2,uni}>1 | fValsPO{3,uni}>1), fValsPO{3,uni}(fValsPO{2,uni}>1 | fValsPO{3,uni}>1));
    poCorr(3,uni) = curPOhtCorr(2,1);
    
    curAllPosF = [fValsPI{1,uni}, fValsPO{1,uni}];
    curAllHeadF = [fValsPI{2,uni}, fValsPO{2,uni}];
    curAllTailF = [fValsPI{3,uni}, fValsPO{3,uni}];
    curALLphCorr = corrcoef(curAllPosF(curAllPosF>1 | curAllHeadF>1), curAllHeadF(curAllPosF>1 | curAllHeadF>1));
    allCorr(1,uni) = curALLphCorr(2,1);
    curALLptCorr = corrcoef(curAllPosF(curAllPosF>1 | curAllTailF>1), curAllTailF(curAllPosF>1 | curAllTailF>1));
    allCorr(2,uni) = curALLptCorr(2,1);
    curALLhtCorr = corrcoef(curAllHeadF(curAllHeadF>1 | curAllTailF>1), curAllTailF(curAllHeadF>1 | curAllTailF>1));
    allCorr(3,uni) = curALLhtCorr(2,1);
    
end
figure; 
sp1 = subplot(3,3,1);
histogram(piCorr(1,:), -1:0.1:1);
hold on;
line(gca, [nanmean(piCorr(1,:)) nanmean(piCorr(1,:))], get(gca, 'ylim'), 'color','r', 'linewidth', 2);
title('PI: Position vs. Head');
sp2 = subplot(3,3,2);
histogram(piCorr(2,:), -1:0.1:1);
hold on;
line(gca, [nanmean(piCorr(2,:)) nanmean(piCorr(2,:))], get(gca, 'ylim'), 'color','r', 'linewidth', 2);
title('PI: Position vs. Tail');
sp3 = subplot(3,3,3);
histogram(piCorr(3,:), -1:0.1:1);
hold on;
line(gca, [nanmean(piCorr(3,:)) nanmean(piCorr(3,:))], get(gca, 'ylim'), 'color','r', 'linewidth', 2);
title('PI: Head vs. Tail');
sp4 = subplot(3,3,4);
histogram(poCorr(1,:), -1:0.1:1);
hold on;
line(gca, [nanmean(poCorr(1,:)) nanmean(poCorr(1,:))], get(gca, 'ylim'), 'color','r', 'linewidth', 2);
title('PO: Position vs. Head');
sp5 = subplot(3,3,5);
histogram(poCorr(2,:), -1:0.1:1);
hold on;
line(gca, [nanmean(poCorr(2,:)) nanmean(poCorr(2,:))], get(gca, 'ylim'), 'color','r', 'linewidth', 2);
title('PO: Position vs. Tail');
sp6 = subplot(3,3,6);
histogram(poCorr(3,:), -1:0.1:1);
hold on;
line(gca, [nanmean(poCorr(3,:)) nanmean(poCorr(3,:))], get(gca, 'ylim'), 'color','r', 'linewidth', 2);
title('PO: Head vs. Tail');
sp7 = subplot(3,3,7);
histogram(allCorr(1,:), -1:0.1:1);
hold on;
line(gca, [nanmean(allCorr(1,:)) nanmean(allCorr(1,:))], get(gca, 'ylim'), 'color','r', 'linewidth', 2);
title('All: Position vs. Head');
sp8 = subplot(3,3,8);
histogram(allCorr(2,:), -1:0.1:1);
hold on;
line(gca, [nanmean(allCorr(2,:)) nanmean(allCorr(2,:))], get(gca, 'ylim'), 'color','r', 'linewidth', 2);
title('All: Position vs. Tail');
sp9 = subplot(3,3,9);
histogram(allCorr(3,:), -1:0.1:1);
hold on;
line(gca, [nanmean(allCorr(3,:)) nanmean(allCorr(3,:))], get(gca, 'ylim'), 'color','r', 'linewidth', 2);
title('All: Head vs. Tail');

linkaxes([sp1, sp2, sp2, sp3, sp4, sp5, sp6, sp7, sp8, sp9], 'y');

orient(gcf, 'tall');
orient(gcf, 'landscape');
annotation('textbox', [0.05 0.9 0.9 0.1], 'String', 'Threshold Fz>1', 'linestyle', 'none', 'FontSize', 20);
annotation('textbox', [0.2 0.9 0.8 0.1], 'String', sprintf('%s Spatial Bin = %.02fmm Gaussian = %.02fms Number of Permutations = %i', seqType, spatialBinSize*1.7, slideWindowSize, numPerms), 'HorizontalAlignment', 'right', 'linestyle', 'none', 'FontSize', 12);
annotation('textbox', 'position', [0.01 0.01 0.9 0.05], 'string', sprintf('%s', cd), 'linestyle', 'none', 'interpreter', 'none');
print('-painters', gcf, '-dpdf', 'LocationVsPositionCorrelationSummary(Thresholded)');
%% Examine the relationship between the peak FR and peak IC
fVals = cellfun(@(a,b)[a,b], fValsPI, fValsPO, 'uniformoutput', 0);
frMean = [mean(pokeInFRisc,3); mean(pokeOutFRisc,3)];
fvMaxNdx = nan(3,length(goodUniNames));
frMaxNdx = nan(1,length(goodUniNames));
for uni = 1:length(goodUniNames)
    fvMaxNdx(1,uni) = find(fVals{1,uni}==max(fVals{1,uni}),1,'first');
    fvMaxNdx(2,uni) = find(fVals{2,uni}==max(fVals{2,uni}),1,'first');
    fvMaxNdx(3,uni) = find(fVals{3,uni}==max(fVals{3,uni}),1,'first');    
    frMaxNdx(uni) = find(frMean(:,uni)==max(frMean(:,uni)),1,'first');
end

figure;
subplot(1,3,1)
corrScatPlot(fvMaxNdx(1,:)', frMaxNdx', 'Index of Max Seq Pos Info', 'Index of Max Firing Rate', []);
subplot(1,3,2)
corrScatPlot(fvMaxNdx(2,:)', frMaxNdx', 'Index of Max Head Loc Info', 'Index of Max Firing Rate', []);
subplot(1,3,3)
corrScatPlot(fvMaxNdx(3,:)', frMaxNdx', 'Index of Max Tail Loc Info', 'Index of Max Firing Rate', []);

orient(gcf, 'tall');
orient(gcf, 'landscape');
annotation('textbox', [0.05 0.9 0.9 0.1], 'String', 'FR vs FZ correlation', 'linestyle', 'none', 'FontSize', 20);
annotation('textbox', [0.2 0.9 0.8 0.1], 'String', sprintf('%s Spatial Bin = %.02fmm Gaussian = %.02fms Number of Permutations = %i', seqType, spatialBinSize*1.7, slideWindowSize, numPerms), 'HorizontalAlignment', 'right', 'linestyle', 'none', 'FontSize', 12);
annotation('textbox', 'position', [0.01 0.01 0.9 0.05], 'string', sprintf('%s', cd), 'linestyle', 'none', 'interpreter', 'none');
print('-painters', gcf, '-dpdf', 'FzVsFRcorr');

figure; 
BarPlotErrorbars([mean(abs(fvMaxNdx(1,:) - frMaxNdx)), mean(abs(fvMaxNdx(2,:) - frMaxNdx)), mean(abs(fvMaxNdx(3,:) - frMaxNdx))],...
    [std(fvMaxNdx(1,:) - frMaxNdx)/sqrt(length(goodUniNames)-1), std(fvMaxNdx(2,:) - frMaxNdx)/sqrt(length(goodUniNames)-1), std(fvMaxNdx(3,:) - frMaxNdx)/sqrt(length(goodUniNames)-1)]);
set(gca, 'xticklabel', {'Position', 'Head', 'Tail'});
title('Temporal Difference From Max Firing Rate');
ylabel('Time (ms)');

orient(gcf, 'tall');
orient(gcf, 'landscape');
annotation('textbox', [0.2 0.9 0.8 0.1], 'String', sprintf('%s Spatial Bin = %.02fmm Gaussian = %.02fms Number of Permutations = %i', seqType, spatialBinSize*1.7, slideWindowSize, numPerms), 'HorizontalAlignment', 'right', 'linestyle', 'none', 'FontSize', 12);
annotation('textbox', 'position', [0.01 0.01 0.9 0.05], 'string', sprintf('%s', cd), 'linestyle', 'none', 'interpreter', 'none');
print('-painters', gcf, '-dpdf', 'PeakInfoFiringRateLatency');
%%
function [fVectRaw, fVectZ] = SlidingFvalCalc(curUniFR, idVect, numPerms)
% UniFR here is organized as a 2-D matrix where each column is a trial and
% each row is a timepoint.
if size(idVect,1) == 1
    idVect = repmat(idVect, [size(curUniFR,1),1]);
end
    
fVectRaw = nan(1,size(curUniFR,1));
parfor t = 1:size(curUniFR,1)
    [~,table,~] = anova1(curUniFR(t,:)', idVect(t,:), 'off');
    if ~isempty(table{2,5})
        fVectRaw(t) = table{2,5};
    else
        fVectRaw(t) = 1;
    end
end
fVectRaw(isnan(fVectRaw)) = 1;
fVectPerm = nan(numPerms, size(curUniFR,1));
for r = 1:numPerms
    idVectShflVect = randperm(size(idVect,2));
    idVectShfl = idVect(:,idVectShflVect);
    parfor t = 1:size(curUniFR,1)
        [~,tableRND,~] = anova1(curUniFR(t,:)', idVectShfl(t,:), 'off');
        if ~isempty(tableRND{2,5})
            fVectPerm(r,t) = tableRND{2,5};
        else
            fVectPerm(r,t) = 1;
        end
    end
end
fVectPerm(isnan(fVectPerm)) = 1;
fVectZ = nan(1,size(curUniFR,1));
for t = 1:size(curUniFR,1)
    tempZ = zscore([fVectRaw(t); fVectPerm(:,t)]);
    fVectZ(t) = tempZ(1);
end

end

%%
function [fVectRaw, fVectZ] = SlidingFvalCalcInteract(curUniFR, idVect1, idVect2, idVect3, numPerms)
if size(idVect1,1) == 1
    idVect1 = repmat(idVect1, [size(curUniFR,1),1]);
end
if size(idVect2,1) == 1
    idVect2 = repmat(idVect2, [size(curUniFR,1),1]);
end
if size(idVect3,1) == 1
    idVect3 = repmat(idVect3, [size(curUniFR,1),1]);
end

fVectRaw = nan(1,size(curUniFR,1));
parfor t = 1:size(curUniFR,1)
    if mean(isnan(idVect3(t,:)))>0.5 || mean(isnan(idVect2(t,:)))>0.5 || mean(isnan(idVect3(t,:)))>0.5
        continue;
    else
        [~,table,~] = anovan(curUniFR(t,:)', {idVect1(t,:), idVect2(t,:), idVect3(t,:)}, 'model', [1 1 1], 'sstype', 'h', 'display', 'off')
        if ~isempty(table{2,5})
            fVectRaw(t) = table{2,6};
        else
            fVectRaw(t) = 1;
        end
    end
end
fVectRaw(isnan(fVectRaw)) = 1;
fVectPerm = nan(numPerms, size(curUniFR,1));
for r = 1:numPerms
    idVect1ShflVect = randperm(size(idVect1,2));
    idVect1Shfl = idVect1(:,idVect1ShflVect);
    idVect2ShflVect = randperm(size(idVect2,2));
    idVect2Shfl = idVect2(:,idVect2ShflVect);
    idVect3ShflVect = randperm(size(idVect3,2));
    idVect3Shfl = idVect1(:,idVect3ShflVect);
    parfor t = 1:size(curUniFR,1)
        if mean(isnan(idVect1Shfl(t,:)))>0.5 || mean(isnan(idVect2Shfl(t,:)))>0.5 || mean(isnan(idVect3Shfl(t,:)))>0.5
        continue;
        else
            [~,table,~] = anovan(curUniFR(t,:)', {idVect1Shfl(t,:), idVect2Shfl(t,:), idVect3Shfl(t,:)}, 'model', [1 1 1], 'sstype', 'h', 'display', 'off')
            if ~isempty(table{2,5})
                fVectPerm(r,t) = table{2,6};
            else
                fVectPerm(r,t) = 1;
            end
        end
    end
end
fVectPerm(isnan(fVectPerm)) = 1;
fVectZ = nan(1,size(curUniFR,1));
for t = 1:size(curUniFR,1)
    tempZ = zscore([fVectRaw(t); fVectPerm(:,t)]);
    fVectZ(t) = tempZ(1);
end

end