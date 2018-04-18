%% F-Ratio Analysis
global spHands
%%
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

%% Analysis Variables
pokeInWindow = [-0.75 0.5];
odorWindow = [-0.75 0.5];
pokeOutWindow = [-0.5 0.75]; 
frontRewardWindow = [-0.5 0.5];
errorWindow = [-0.5 0.5];

pehBinSize = 0.125;
fMax = 5;

%% Load Relevant Data
load(fileNames{cellfun(@(a)~isempty(a), strfind(fileNames, 'BehaviorMatrix'))});
load(fileNames{cellfun(@(a)~isempty(a), strfind(fileNames, 'EnsembleMatrix'))});

%% Select Units
spr = [ensembleUnitSummaries.Spike_Phase_Relations];
rip = [spr.Ripple];
ripR = [rip.R_Length];
uniLog = ripR>0.25;

curUniMtx = ensembleMatrix(:,2:end);
curUniMtx = curUniMtx(:,uniLog);
%% Calculate Event F-Ratios
figure;
% Poke In
pokeInBehavMatrix = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, pokeInWindow, 'PokeIn');
perfLog = [pokeInBehavMatrix.Performance] == 1;
firstPosLog = [pokeInBehavMatrix.Position] == 1;
pokeInFvalsCorr = CalculatePosVsItmFrats_SM(curUniMtx, pokeInBehavMatrix(perfLog&~firstPosLog), pokeInWindow, pehBinSize, fMax);
pokeInFvalsInCorr = CalculatePosVsItmFrats_SM(curUniMtx, pokeInBehavMatrix(~perfLog&~firstPosLog), pokeInWindow, pehBinSize, fMax);

pokeInXTicks = (pokeInWindow(1):pehBinSize:(pokeInWindow(2)-pehBinSize))+(pehBinSize/2);
PlotFvalsFig(pokeInFvalsCorr, pokeInFvalsInCorr, 1:3, pokeInXTicks, 'Poke In');
drawnow

% Odor
odorDeliveryBehavMatrix = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, odorWindow, 'Odor');
perfLog = [odorDeliveryBehavMatrix.Performance] == 1;
firstPosLog = [odorDeliveryBehavMatrix.Position] == 1;
odorFvalsCorr = CalculatePosVsItmFrats_SM(curUniMtx, odorDeliveryBehavMatrix(perfLog&~firstPosLog), odorWindow, pehBinSize, fMax);
odorFvalsInCorr = CalculatePosVsItmFrats_SM(curUniMtx, odorDeliveryBehavMatrix(~perfLog&~firstPosLog), odorWindow, pehBinSize, fMax);

odorXTicks = (odorWindow(1):pehBinSize:(odorWindow(2)-pehBinSize))+(pehBinSize/2);
PlotFvalsFig(odorFvalsCorr, odorFvalsInCorr, 4:6, odorXTicks, 'Odor Delivery');
drawnow

% PokeOut
pokeOutBehavMatrix = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, pokeOutWindow, 'PokeOut');
perfLog = [pokeOutBehavMatrix.Performance] == 1;
firstPosLog = [pokeOutBehavMatrix.Position] == 1;
pokeOutFvalsCorr = CalculatePosVsItmFrats_SM(curUniMtx, pokeOutBehavMatrix(perfLog&~firstPosLog), pokeOutWindow, pehBinSize, fMax);
pokeOutFvalsInCorr = CalculatePosVsItmFrats_SM(curUniMtx, pokeOutBehavMatrix(~perfLog&~firstPosLog), pokeOutWindow, pehBinSize, fMax);

pokeOutXTicks = (pokeOutWindow(1):pehBinSize:(pokeOutWindow(2)-pehBinSize))+(pehBinSize/2);
PlotFvalsFig(pokeOutFvalsCorr, pokeOutFvalsInCorr, 7:9, pokeOutXTicks, 'Poke Withdrawal');
drawnow

% FrontReward
frontRewardBehavMatrix = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, frontRewardWindow, 'FrontReward');
perfLog = [frontRewardBehavMatrix.Performance] == 1;
firstPosLog = [frontRewardBehavMatrix.Position] == 1;
frontRewardFvalsCorr = CalculatePosVsItmFrats_SM(curUniMtx, frontRewardBehavMatrix(perfLog&~firstPosLog), frontRewardWindow, pehBinSize, fMax);
frontRewardFvalsInCorr = CalculatePosVsItmFrats_SM(curUniMtx, frontRewardBehavMatrix(~perfLog&~firstPosLog), frontRewardWindow, pehBinSize, fMax);

frontRewardXTicks = (frontRewardWindow(1):pehBinSize:(frontRewardWindow(2)-pehBinSize))+(pehBinSize/2);
PlotFvalsFig(frontRewardFvalsCorr, frontRewardFvalsInCorr, 10:12, frontRewardXTicks, 'Front Reward Delivery');
drawnow

% Error
errorBehavMatrix = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, errorWindow, 'ErrorSignal');
perfLog = [errorBehavMatrix.Performance] == 1;
firstPosLog = [errorBehavMatrix.Position] == 1;
errorFvalsCorr = CalculatePosVsItmFrats_SM(curUniMtx, errorBehavMatrix(perfLog&~firstPosLog), errorWindow, pehBinSize, fMax);
errorFvalsInCorr = CalculatePosVsItmFrats_SM(curUniMtx, errorBehavMatrix(~perfLog&~firstPosLog), errorWindow, pehBinSize, fMax);

errorXTicks = (errorWindow(1):pehBinSize:(errorWindow(2)-pehBinSize))+(pehBinSize/2);
PlotFvalsFig(errorFvalsCorr, errorFvalsInCorr, 13:15, errorXTicks, 'Error Signal');
drawnow

linkaxes([spHands(1:3:13), spHands(2:3:14)], 'y');
linkaxes(spHands(3:3:15), 'y');

for spNum = 1:15
    subplot(5,3,spNum)
    if mod(spNum,3)==0
        line(get(gca, 'xlim'), [0 0], 'color', 'black', 'linestyle', ':', 'linewidth', 1.5);
        line([0 0], get(gca, 'ylim'), 'color', 'black', 'linestyle', ':', 'linewidth', 1.5);
    else
        line(get(gca, 'xlim'), [1 1], 'color', 'black', 'linestyle', ':', 'linewidth', 1.5);
        line([0 0], get(gca, 'ylim'), 'color', 'black', 'linestyle', ':', 'linewidth', 1.5);
    end
end
%%
function PlotFvalsFig(fValMatrixCorr, fValMatrixInCorr, fValSubplots, xTicks, eventID)
global spHands
    spHands(fValSubplots(1)) = subplot(5,3,fValSubplots(1));
    hands1 = PlotFvals(fValMatrixCorr, 1, xTicks, 'Red');
    hold on;
    hands2 = PlotFvals(fValMatrixCorr, 2, xTicks, 'Blue');
    axis tight;
    legend([hands1 hands2], 'Position', 'Odor');
    title([eventID, ' Correct']);
    xlabel('Time Relative to Event');
    ylabel('F-Ratio');
    
    spHands(fValSubplots(2)) = subplot(5,3,fValSubplots(2));
    hands3 = PlotFvals(fValMatrixInCorr, 1, xTicks, 'Red');
    hold on;
    hands4 = PlotFvals(fValMatrixInCorr, 2, xTicks, 'Blue');
    axis tight;
    legend([hands3 hands4], 'Position', 'Odor');
    title([eventID, ' Incorrect']);
    xlabel('Time Relative to Event');
    ylabel('F-Ratio');
    
    spHands(fValSubplots(3)) = subplot(5,3,fValSubplots(3));
    hands5 = PlotFvals(fValMatrixCorr, 3, xTicks, 'Green');
    hold on;
    hands6 = PlotFvals(fValMatrixInCorr, 3, xTicks, 'Black');
    axis tight;
    legend([hands5 hands6], 'Correct', 'Incorrect');
    title([{eventID}; 'Correct vs Incorrect']);
    xlabel('Time Relative to Event');
    ylabel({'F-Ratio'; 'Difference'});
end

function axHand = PlotFvals(fValMatrix, fValColID, xTicks, color)
    fVals = cell2mat(fValMatrix(:,fValColID));
    if ~isempty(fVals)
        axHand = PlotLineAndFilledError(xTicks, nanmean(fVals),...
            nanstd(fVals,0,1)./(sum(~isnan(fVals))-1), color);
    else
        axHand = [];
    end
end