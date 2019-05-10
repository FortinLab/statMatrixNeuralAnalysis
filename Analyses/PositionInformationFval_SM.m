
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
slideWindowSize = 30;
spatialBinSize = 5;
numPerms = 50;

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
uniInstFR = [ensembleMatrix(:,1), nan(size(ensembleMatrix,1), size(ensembleMatrix,2)-1)];
for uni = 2:size(ensembleMatrix,2)
    uniInstFR(:,uni) = conv(ensembleMatrix(:,uni), instFRgauss, 'same');
end
% 
% % Alternatively, just bin the firing rates
% uniInstFR = [ensembleMatrix(:,1), nan(size(ensembleMatrix,1), size(ensembleMatrix,2)-1)];
% for uni = 2:size(ensembleMatrix,2)
%     uniInstFR(:,uni) = conv(ensembleMatrix(:,uni), ones(1,slideWindowSize)/(length(instFRgauss)*mode(diff(behavMatrix(:,1)))), 'same'); 
% end

%% Organize Trial Data
% Create Extraction Matrices
pokeInTrialMatrix = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [-0.2 0.5], 'PokeIn');
pokeInTSs = behavMatrix(pokeInTrialMatrix(1).TrialLogVect,1)-behavMatrix(pokeInTrialMatrix(1).PokeInIndex,1);
pokeOutTrialMatrix = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [-0.5 0.2], 'PokeOut');
pokeOutTSs = behavMatrix(pokeOutTrialMatrix(1).TrialLogVect,1)-behavMatrix(pokeOutTrialMatrix(1).PokeOutIndex,1);

% Create Trial Logical Vectors
corrTrlLog = [pokeInTrialMatrix.Performance];
isLog = [pokeInTrialMatrix.TranspositionDistance]==0;
posIDsISC = [pokeInTrialMatrix(corrTrlLog & isLog).Position];

% Pull out spiking activity
pokeInFR = ExtractTrialData_SM(pokeInTrialMatrix, uniInstFR); %#ok<*NODEF>
pokeInFRisc = cell2mat(reshape(pokeInFR(corrTrlLog & isLog), [1,1,sum(corrTrlLog & isLog)]));
pokeOutFR = ExtractTrialData_SM(pokeOutTrialMatrix, uniInstFR);
pokeOutFRisc = cell2mat(reshape(pokeOutFR(corrTrlLog & isLog), [1,1,sum(corrTrlLog & isLog)]));

% Pull out and interpolate positional information
pokeInOrientation = ExtractTrialData_SM(pokeInTrialMatrix, orientMatrix); %#ok<*NODEF>
pokeInOrientationInterp = cell(size(pokeInOrientation));
pokeOutOrientation = ExtractTrialData_SM(pokeOutTrialMatrix, orientMatrix);
pokeOutOrientationInterp = cell(size(pokeOutOrientation));

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
pokeInOrientISC = cell2mat(reshape(pokeInOrientationInterp(corrTrlLog & isLog), [1,1,sum(corrTrlLog & isLog)]));
pokeOutOrientISC = cell2mat(reshape(pokeOutOrientationInterp(corrTrlLog & isLog), [1,1,sum(corrTrlLog & isLog)]));

pokeInOrientBinsISC = cell2mat(reshape(orientBinValsPI(corrTrlLog & isLog), [1,1,sum(corrTrlLog & isLog)]));
pokeOutOrientBinsISC = cell2mat(reshape(orientBinValsPO(corrTrlLog & isLog), [1,1,sum(corrTrlLog & isLog)]));

%% Run Sliding F-Ratio Analysis
for uni = 1:length(ensembleUnitSummaries)
    curUniPokeInFR = reshape(pokeInFRisc(:,uni+1,:), [size(pokeInFRisc,1), size(pokeInFRisc,3), 1]);
    curUniPokeOutFR = reshape(pokeOutFRisc(:,uni+1,:), [size(pokeOutFRisc,1), size(pokeOutFRisc,3), 1]);
    
    % Poke In Oriented
    [curUniPIfvalPOS, curUniPIfvalPOSz] = SlidingFvalCalc(curUniPokeInFR, posIDsISC, numPerms);
    [curUniPIfvalHeadLoc, curUniPIfvalHeadLocZ] = SlidingFvalCalc(curUniPokeInFR, reshape(pokeInOrientBinsISC(:,1,:), [size(pokeInOrientBinsISC,1), size(pokeInOrientBinsISC,3), 1]), numPerms);
    [curUniPIfvalTailLoc, curUniPIfvalTailLocZ] = SlidingFvalCalc(curUniPokeInFR, reshape(pokeOutOrientBinsISC(:,2,:), [size(pokeOutOrientBinsISC,1), size(pokeOutOrientBinsISC,3), 1]), numPerms);

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
    [curUniPOfvalHeadLoc, curUniPOfvalHeadLocZ] = SlidingFvalCalc(curUniPokeInFR, reshape(pokeOutOrientBinsISC(:,1,:), [size(pokeOutOrientBinsISC,1), size(pokeOutOrientBinsISC,3), 1]), numPerms);
    [curUniPOfvalTailLoc, curUniPOfvalTailLocZ] = SlidingFvalCalc(curUniPokeInFR, reshape(pokeOutOrientBinsISC(:,2,:), [size(pokeOutOrientBinsISC,1), size(pokeOutOrientBinsISC,3), 1]), numPerms);

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
    legend('Sequence Position', 'Head', 'Tail', 'location', 'best');
%     legend('Sequence Position', 'Head Loc X', 'Head Loc Y', 'Tail Loc X', 'Tail Loc Y', 'Port Angle', 'Head Angle', 'Tail Angle', 'Head-Tail', 'Head-Port', 'Tail-Port', 'location', 'best');
    
    sp2 = subplot(2,1,2);
    plot(pokeOutTSs, curUniPOfvalPOSz, 'linewidth', 2);
    hold on;
    plot(pokeOutTSs, curUniPOfvalHeadLocZ, 'k');
    plot(pokeOutTSs, curUniPOfvalTailLocZ, '--k');
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
    annotation('textbox', [0.05 0.9 0.9 0.1], 'String', ensembleMatrixColIDs{uni+1}, 'linestyle', 'none', 'FontSize', 20);
    annotation('textbox', [0.2 0.9 0.8 0.1], 'String', sprintf('Spatial Bin = %.02fmm Gaussian = %.02fms Number of Permutations = %i', spatialBinSize*1.7, slideWindowSize, numPerms), 'HorizontalAlignment', 'right', 'linestyle', 'none', 'FontSize', 12);
    annotation('textbox', 'position', [0.01 0.01 0.9 0.05], 'string', sprintf('%s', cd), 'linestyle', 'none', 'interpreter', 'none');
    
    orient(gcf, 'tall');
    orient(gcf, 'landscape');
%     print
    print('-painters', gcf, '-dpdf', sprintf('%s_OrientFvalInfo_Summary', ensembleMatrixColIDs{uni+1}));        
end

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