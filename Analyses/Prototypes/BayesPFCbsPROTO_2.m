%%
clc
clear all

%% Runtime variables
binSize = 20;
dsRate = 5;
pokeInWindows = [0 1.2];
pokeOutWindows = [-1.2 0];

%%
smPath = uigetdir;
cd(smPath);
files = dir(smPath);
fileNames = {files.name};
% Load the behavior matrix file for poke events plots
behMatFile = fileNames{cellfun(@(a)~isempty(a), strfind(fileNames, 'BehaviorMatrix'))};
nsmblMatFile = fileNames{cellfun(@(a)~isempty(a), strfind(fileNames, 'EnsembleMatrix'))};
load([smPath '\' behMatFile]);
load([smPath '\' nsmblMatFile]);
% Identify list of all statMatrix files
smFileList = fileNames(cellfun(@(a)~isempty(a), regexp(fileNames, '_SM\>')))';

%% Extract Behavioral Periods
% Taking 1/2 the binSize on either end to get rid of edge effects.
pokeInTrialPeriodTD = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [pokeInWindows(1)-(binSize/2/1000) pokeInWindows(2)+(binSize/2/1000)], 'PokeIn');
pokeInTimes = behavMatrix(pokeInTrialPeriodTD(1).TrialLogVect,1) - behavMatrix(pokeInTrialPeriodTD(1).PokeInIndex,1);
pokeOutTrialPeriodTD = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [pokeOutWindows(1)-(binSize/2/1000) pokeOutWindows(2)+(binSize/2/1000)], 'PokeOut');
pokeOutTimes = behavMatrix(pokeOutTrialPeriodTD(1).TrialLogVect,1) - behavMatrix(pokeOutTrialPeriodTD(1).PokeOutIndex,1);

% pokeInTrialPeriodTD = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [pokeInWindows(1)-(binSize/1000) pokeInWindows(2)+(binSize/1000)], 'PokeIn');
% pokeInTimes = behavMatrix(pokeInTrialPeriodTD(1).TrialLogVect,1) - behavMatrix(pokeInTrialPeriodTD(1).PokeInIndex,1);
% pokeOutTrialPeriodTD = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [pokeOutWindows(1)-(binSize/1000) pokeOutWindows(2)+(binSize/1000)], 'PokeOut');
% pokeOutTimes = behavMatrix(pokeOutTrialPeriodTD(1).TrialLogVect,1) - behavMatrix(pokeOutTrialPeriodTD(1).PokeOutIndex,1);

pokeInEnsemble = ExtractTrialData_SM(pokeInTrialPeriodTD, ensembleMatrix(:,2:end)); %#ok<*NODEF>
pokeInEnsembleMtx = cell2mat(reshape(pokeInEnsemble, [1 1 length(pokeInEnsemble)]));
pokeOutEnsemble = ExtractTrialData_SM(pokeOutTrialPeriodTD, ensembleMatrix(:,2:end)); %#ok<*NODEF>
pokeOutEnsembleMtx = cell2mat(reshape(pokeOutEnsemble, [1 1 length(pokeOutEnsemble)]));


%% Bin the spiking data
% First convolve the entire trialEnsembleMtx with a square to bin the
% spikes
pokeInBinnedEnsembleMtx = nan(size(pokeInEnsembleMtx));
pokeOutBinnedEnsembleMtx = nan(size(pokeOutEnsembleMtx));
for t = 1:size(pokeInEnsembleMtx,3)
    for u = 1:size(pokeInEnsembleMtx,2)
        pokeInBinnedEnsembleMtx(:,u,t) = conv(pokeInEnsembleMtx(:,u,t), ones(1,binSize)./(binSize/1000), 'same');
        pokeOutBinnedEnsembleMtx(:,u,t) = conv(pokeOutEnsembleMtx(:,u,t), ones(1,binSize)./(binSize/1000), 'same');
    end
end
% Now remove the binSize/2 padding
unPaddedPokeInBinnedEnsembleMtx = pokeInBinnedEnsembleMtx((binSize/2)+1:end-(binSize/2),:,:);
pokeInTimesUnpadded = pokeInTimes((binSize/2)+1:end-(binSize/2));
unPaddedPokeOutBinnedEnsembleMtx = pokeOutBinnedEnsembleMtx((binSize/2)+1:end-(binSize/2),:,:);
pokeOutTimesUnpadded = pokeOutTimes((binSize/2)+1:end-(binSize/2));
% Now downsample the binned matrix
dsVect = downsample(1:size(unPaddedPokeInBinnedEnsembleMtx,1), dsRate);
pokeInSpikeMatrix = unPaddedPokeInBinnedEnsembleMtx(dsVect,:,:);
pokeOutSpikeMatrix = unPaddedPokeOutBinnedEnsembleMtx(dsVect,:,:);
trialTimePI = pokeInTimesUnpadded(dsVect);
trialTimePO = pokeOutTimesUnpadded(dsVect);

%% Alternative way to bin the spike data
% mansWin1 = 1:dsRate:size(pokeInEnsembleMtx,1)-binSize;
% mansWin2 = mansWin1+binSize;
% 
% pokeInBinnedEnsembleMtx = nan(length(mansWin1), size(pokeInEnsembleMtx,2), size(pokeInEnsembleMtx,3));
% pokeOutBinnedEnsembleMtx = nan(length(mansWin1), size(pokeOutEnsembleMtx,2), size(pokeOutEnsembleMtx,3));
% for trl = 1:size(pokeInEnsembleMtx,3)
%     for u = 1:size(pokeInEnsembleMtx,2)        
%         curUniPI = pokeInEnsembleMtx(:,u,trl);
%         pokeInBinnedEnsembleMtx(:,u,trl) = cell2mat(arrayfun(@(a,b)sum(curUniPI(a:b)), mansWin1, mansWin2, 'uniformoutput', 0));
%         curUniPO = pokeOutEnsembleMtx(:,u,trl);
%         pokeOutBinnedEnsembleMtx(:,u,trl) = cell2mat(arrayfun(@(a,b)sum(curUniPO(a:b)), mansWin1, mansWin2, 'uniformoutput', 0));
%     end
% end
% pokeInBinnedEnsembleMtx = pokeInBinnedEnsembleMtx./(binSize/1000);
% pokeOutBinnedEnsembleMtx = pokeOutBinnedEnsembleMtx./(binSize/1000);
% pokeInTimesDS = pokeInTimes(mansWin1);
% pokeInUnPadLog = pokeInTimesDS>=pokeInWindows(1) & pokeInTimesDS<=pokeInWindows(2);
% trialTimePI = pokeInTimesDS(pokeInUnPadLog);
% pokeInSpikeMatrix = pokeInBinnedEnsembleMtx(pokeInUnPadLog,:,:);
% pokeOutTimesDS = pokeOutTimes(mansWin1);
% pokeOutUnPadLog = pokeOutTimesDS>=pokeOutWindows(1) & pokeOutTimesDS<=pokeOutWindows(2);
% trialTimePO = pokeOutTimesDS(pokeOutUnPadLog);
% pokeOutSpikeMatrix = pokeOutBinnedEnsembleMtx(pokeOutUnPadLog,:,:);


%% Create Logical Vectors
perfLog = [pokeInTrialPeriodTD.Performance];
inSeqLog = [pokeInTrialPeriodTD.TranspositionDistance]==0;
outSeqLog = [pokeInTrialPeriodTD.TranspositionDistance]~=0 & abs([pokeInTrialPeriodTD.TranspositionDistance])<10;
odorAlog = [pokeInTrialPeriodTD.Odor] == 1;
odorBlog = [pokeInTrialPeriodTD.Odor] == 2;
odorClog = [pokeInTrialPeriodTD.Odor] == 3;
odorDlog = [pokeInTrialPeriodTD.Odor] == 4;

fullInSeqSeqsStart = find(conv([pokeInTrialPeriodTD.Odor], 1:4, 'valid')==20 & conv([pokeInTrialPeriodTD.Position], 1:4, 'valid')==20 & conv([pokeInTrialPeriodTD.Performance]*1, ones(1,4), 'valid')==4);
inSeqSeqs = nan(3,length(fullInSeqSeqsStart));
for iS = 1:length(fullInSeqSeqsStart)
    inSeqSeqs(1,iS) = fullInSeqSeqsStart(iS);
    inSeqSeqs(2,iS) = fullInSeqSeqsStart(iS) + 1;
    inSeqSeqs(3,iS) = fullInSeqSeqsStart(iS) + 2;
    inSeqSeqs(4,iS) = fullInSeqSeqsStart(iS) + 3;
end
fullInSeqLog = false(1,length(pokeInTrialPeriodTD));
fullInSeqLog(inSeqSeqs(:)) = true;

%% 
uniFRthreshLog = max([mean(pokeInSpikeMatrix,3); mean(pokeOutSpikeMatrix,3)])<1;
pokeInSpkMtx = pokeInSpikeMatrix;
pokeInSpkMtx(:,uniFRthreshLog,:) = [];
pokeOutSpkMtx = pokeOutSpikeMatrix;
pokeOutSpkMtx(:,uniFRthreshLog,:) = [];
goodUniNames = {ensembleUnitSummaries(~uniFRthreshLog).UnitName};

%%
figure
corrISmtxPI = mean(pokeInSpkMtx(:,:,perfLog & fullInSeqLog),3);                         % All InSeq Correct Trials
[postPI] = CalculatePostProb(corrISmtxPI, pokeInSpkMtx(:,:,perfLog & inSeqLog & ~fullInSeqLog), binSize);
subplot(4,4,[1,2,5,6])
PlotPostMtx(trialTimePI, postPI, 'InSeq Correct Trials (Poke In)');
set(gca, 'clim', [0 0.5]);

corrISmtxPO = mean(pokeOutSpkMtx(:,:,perfLog & fullInSeqLog),3);
[postPO] = CalculatePostProb(corrISmtxPO, pokeOutSpkMtx(:,:,perfLog & inSeqLog & ~fullInSeqLog), binSize);
subplot(4,4,[9,10,13,14])
PlotPostMtx(trialTimePO, postPO, 'InSeq Correct Trials (Poke Out)');
set(gca, 'clim', [0 0.5]);


nonAIStrialData = pokeInTrialPeriodTD(perfLog & inSeqLog & ~fullInSeqLog);
nonAISodors = [nonAIStrialData.Odor];
nonAIStrialNums = [nonAIStrialData.TrialNum];

preOStrialLog = false(size(nonAIStrialNums));
postOStrialLog = false(size(nonAIStrialNums));
for trl = 1:length(nonAIStrialNums)
    if nonAIStrialNums(trl)~=1 && nonAIStrialNums(trl) ~= length(pokeInTrialPeriodTD)
        prevTrial = pokeInTrialPeriodTD(nonAIStrialNums(trl)-1);
        curTrl = pokeInTrialPeriodTD(nonAIStrialNums(trl));
        nextTrial = pokeInTrialPeriodTD(nonAIStrialNums(trl)+1);
        if curTrl.Position ~= 1
            if curTrl.Position - prevTrial.Position == 1 && prevTrial.TranspositionDistance ~= 0
                postOStrialLog(trl) = true;
                preOStrialLog(trl-2) = true;
            end
        end
   end
end

preOSpi = postPI(:,:,preOStrialLog);
subplot(4,4,3)
PlotPostMtx(trialTimePI, preOSpi, 'Before OutSeq');
set(gca, 'clim', [0 0.5]);
postOSpi = postPI(:,:,postOStrialLog);
subplot(4,4,7)
PlotPostMtx(trialTimePI, postOSpi, 'After OutSeq');
set(gca, 'clim', [0 0.5]);
postPreDiffPI = postOSpi - preOSpi;
subplot(4,4,4)
PlotPostMtx(trialTimePI, postPreDiffPI, 'After-Before');
set(gca, 'clim', [-0.25 0.25]);


preOSpo = postPO(:,:,preOStrialLog);
subplot(4,4,11)
PlotPostMtx(trialTimePO, preOSpo, 'Before OutSeq');
set(gca, 'clim', [0 0.5]);
postOSpo = postPO(:,:,postOStrialLog);
subplot(4,4,15)
PlotPostMtx(trialTimePO, postOSpo, 'After OutSeq');
set(gca, 'clim', [0 0.5]);
postPreDiffPO = postOSpo - preOSpo;
subplot(4,4,12)
PlotPostMtx(trialTimePO, postPreDiffPO, 'After-Before');
set(gca, 'clim', [-0.25 0.25]);

annotation('textbox', 'position', [0.5 0.935 0.5 0.05], 'String', ['\bf\fontsize{10}' sprintf('Bin = %i ms; Step = %i ms', binSize, dsRate)],...
    'linestyle', 'none', 'horizontalalignment', 'right');
curDir = cd;
annotation('textbox', 'position', [0.025 0.025 0.7 0.05], 'String', curDir,...
    'linestyle', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
colormap jet
axesHandles = findobj(get(gcf,'Children'), 'flat','Type','axes');
axis(axesHandles,'square')
orient(gcf, 'tall');
orient(gcf, 'landscape');
drawnow;




%%
function cAx = PlotPostMtx(trialTimes, postMtx, id)
imagesc(trialTimes, trialTimes, nanmean(postMtx,3));
set(gca, 'ydir', 'normal')
hold on;
% line(trialTimes, trialTimes, 'linestyle', ':', 'color', 'w', 'linewidth', 2);
title(id);
xlabel('True Time');
ylabel('Decoded Time');
drawnow
cAx = get(gca, 'clim');
end

%%
function [postNorm, postRaw] = CalculatePostProb(prior, obsv, binSize)
propVect = 1./sum(prior,2);
postNorm = nan(size(obsv,1), size(obsv,1), size(obsv,3));
postRaw = nan(size(obsv,1), size(obsv,1), size(obsv,3));
for trl = 1:size(obsv,3)
    for t = 1:size(prior,1)
        p = nan(size(prior));
        curPopVect = obsv(t,:,trl);
        curPopFact = factorial(curPopVect);
        for u = 1:size(prior,2)
            curAvgUniFR = prior(:,u);
            p(:,u) = (((binSize/1000)*curAvgUniFR).^curPopVect(u))./curPopFact(u);
        end        
        pp = prod(p,2);
        ee = exp(-(binSize/1000*sum(prior,2)));
        tempPost = propVect.*pp.*ee;
        postRaw(t,:,trl) = tempPost;
%         postNorm(t,:,trl) = tempPost./max(tempPost);                       % Probably wrong
        postNorm(t,:,trl) = tempPost==max(tempPost);
    end
end
end