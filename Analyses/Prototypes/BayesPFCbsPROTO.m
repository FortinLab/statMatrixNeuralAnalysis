%%
clc
clear all

%% Runtime variables
binSize = 100;
dsRate = 5;

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
trialPeriodTD = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [-0.5-(binSize/2/1000) 1.5+(binSize/2/1000)], 'PokeIn');
% trialPeriodTD = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [-1.5-(binSize/2/1000) 0.5+(binSize/2/1000)], 'PokeOut');
trialEnsemble = ExtractTrialData_SM(trialPeriodTD, ensembleMatrix(:,2:end)); %#ok<*NODEF>
trialEnsembleMtx = cell2mat(reshape(trialEnsemble, [1 1 length(trialEnsemble)]));

trialTimes = behavMatrix(trialPeriodTD(1).TrialLogVect,1) - behavMatrix(trialPeriodTD(1).PokeInIndex,1);
% trialTimes = behavMatrix(trialPeriodTD(1).TrialLogVect,1) - behavMatrix(trialPeriodTD(1).PokeOutIndex,1);

clear behavMatrix ensembleMatrix
%% Bin the spiking data
% First convolve the entire trialEnsembleMtx with a square to bin the
% spikes
binnedEnsembleMtx = nan(size(trialEnsembleMtx));
for t = 1:size(trialEnsembleMtx,3)
    for u = 1:size(trialEnsembleMtx,2)
        binnedEnsembleMtx(:,u,t) = conv(trialEnsembleMtx(:,u,t), ones(1,binSize)./(binSize/1000), 'same');
%         binnedEnsembleMtx(:,u,t) = conv(trialEnsembleMtx(:,u,t), ones(1,binSize), 'same');
    end
end
% Now remove the binSize/2 padding
unPaddedBinnedEnsembleMtx = binnedEnsembleMtx((binSize/2)+1:end-(binSize/2),:,:);
trialTimes = trialTimes((binSize/2)+1:end-(binSize/2));
% Now downsample the binned matrix
dsVect = downsample(1:size(unPaddedBinnedEnsembleMtx,1), dsRate);
spikeMatrix = unPaddedBinnedEnsembleMtx(dsVect,:,:);
trialTime = trialTimes(dsVect);

clear trialEnsembleMtx binnedEnsembleMtx unPaddedBinnedEnsembleMtx
%% Create Logical Vectors
perfLog = [trialPeriodTD.Performance];
inSeqLog = [trialPeriodTD.TranspositionDistance]==0;
outSeqLog = [trialPeriodTD.TranspositionDistance]~=0 & abs([trialPeriodTD.TranspositionDistance])<10;
odorAlog = [trialPeriodTD.Odor] == 1;
odorBlog = [trialPeriodTD.Odor] == 2;
odorClog = [trialPeriodTD.Odor] == 3;
odorDlog = [trialPeriodTD.Odor] == 4;

fullInSeqSeqsStart = find(conv([trialPeriodTD.Odor], 1:4, 'valid')==20 & conv([trialPeriodTD.Position], 1:4, 'valid')==20 & conv([trialPeriodTD.Performance]*1, ones(1,4), 'valid')==4);
inSeqSeqs = nan(3,length(fullInSeqSeqsStart));
for iS = 1:length(fullInSeqSeqsStart)
    inSeqSeqs(1,iS) = fullInSeqSeqsStart(iS);
    inSeqSeqs(2,iS) = fullInSeqSeqsStart(iS) + 1;
    inSeqSeqs(3,iS) = fullInSeqSeqsStart(iS) + 2;
    inSeqSeqs(4,iS) = fullInSeqSeqsStart(iS) + 3;
end
fullInSeqLog = false(1,length(trialPeriodTD));
fullInSeqLog(inSeqSeqs(:)) = true;

nonAIStrialData = trialPeriodTD(perfLog & inSeqLog & ~fullInSeqLog);
%% 
uniFRthreshLog = max(mean(spikeMatrix,3))<1;
spkMtx = spikeMatrix;
spkMtx(:,uniFRthreshLog,:) = [];
goodUniNames = {ensembleUnitSummaries(~uniFRthreshLog).UnitName};

clear spikeMatrix
%% Decode Trial Time Across Odors
figure;
corrAisMtx = mean(spkMtx(:,:,perfLog & fullInSeqLog & odorAlog),3);             % All A InSeq Correct Trials
[postAnorm, postAraw] = CalculatePostProb(corrAisMtx, spkMtx(:,:,perfLog & inSeqLog & odorAlog & ~fullInSeqLog), binSize);
subplot(2,2,1);
aCaxis = PlotPostMtx(trialTimes, postAnorm, 'Odor A');

corrBisMtx = mean(spkMtx(:,:,perfLog & fullInSeqLog & odorBlog),3);             % All B InSeq Correct Trials
[postBnorm, postBraw] = CalculatePostProb(corrBisMtx, spkMtx(:,:,perfLog & inSeqLog & odorBlog & ~fullInSeqLog), binSize);
subplot(2,2,2);
bCaxis = PlotPostMtx(trialTimes, postBnorm, 'Odor B');

corrCisMtx = mean(spkMtx(:,:,perfLog & fullInSeqLog & odorClog),3);             % All C InSeq Correct Trials
[postCnorm, postCraw] = CalculatePostProb(corrCisMtx, spkMtx(:,:,perfLog & inSeqLog & odorClog & ~fullInSeqLog), binSize);
subplot(2,2,3);
cCaxis = PlotPostMtx(trialTimes, postCnorm, 'Odor C');

corrDisMtx = mean(spkMtx(:,:,perfLog & fullInSeqLog & odorDlog),3);             % All D InSeq Correct Trials
[postDnorm, postDraw] = CalculatePostProb(corrDisMtx, spkMtx(:,:,perfLog & inSeqLog & odorDlog & ~fullInSeqLog), binSize);
subplot(2,2,4);
dCaxis = PlotPostMtx(trialTimes, postDnorm, 'Odor D');

cAx = [min([aCaxis, bCaxis, cCaxis, dCaxis]), max([aCaxis, bCaxis, cCaxis, dCaxis])/3];

annotation('textbox', 'position', [0.5 0.935 0.5 0.05], 'String', ['\bf\fontsize{10}' sprintf('Bin = %i ms; Step = %i ms', binSize, dsRate)],...
    'linestyle', 'none', 'horizontalalignment', 'right');
curDir = cd;
annotation('textbox', 'position', [0.025 0.025 0.7 0.05], 'String', curDir,...
    'linestyle', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
colormap jet
axesHandles = findobj(get(gcf,'Children'), 'flat','Type','axes');
axis(axesHandles,'square')
set(axesHandles, 'clim', cAx);
orient(gcf, 'tall');
orient(gcf, 'landscape');

%% Decode Trial Time and Odor Across Odors
nonAIStrials = spkMtx(:,:,perfLog & inSeqLog & ~fullInSeqLog);
nonAISodors = [nonAIStrialData.Odor];
corrAisMtx = mean(spkMtx(:,:,perfLog & fullInSeqLog & odorAlog),3);             % All A InSeq Correct Trials Prior
corrBisMtx = mean(spkMtx(:,:,perfLog & fullInSeqLog & odorBlog),3);             % All B InSeq Correct Trials Prior
corrCisMtx = mean(spkMtx(:,:,perfLog & fullInSeqLog & odorClog),3);             % All C InSeq Correct Trials Prior
corrDisMtx = mean(spkMtx(:,:,perfLog & fullInSeqLog & odorDlog),3);             % All D InSeq Correct Trials Prior

figure;
oAfr = subplot(4,4,1);
imagesc(trialTimes, 1:size(spkMtx,2), corrAisMtx');
title('Odor A');
abCorr = subplot(4,4,2);
corrScatPlot(corrAisMtx(:),corrBisMtx(:), 'A','B', 'markerStyle', '.', 'markerColor', 'k');
acCorr = subplot(4,4,3);
corrScatPlot(corrAisMtx(:),corrCisMtx(:), 'A','C', 'markerStyle', '.', 'markerColor', 'k');
adCorr = subplot(4,4,4);
corrScatPlot(corrAisMtx(:),corrDisMtx(:), 'A','B', 'markerStyle', '.', 'markerColor', 'k');

oBfr = subplot(4,4,6);
imagesc(trialTimes, 1:size(spkMtx,2), corrBisMtx');
title('Odor B');
bcCorr = subplot(4,4,7);
corrScatPlot(corrBisMtx(:),corrCisMtx(:), 'B','C', 'markerStyle', '.', 'markerColor', 'k');
bdCorr = subplot(4,4,8);
corrScatPlot(corrBisMtx(:),corrDisMtx(:), 'B','D', 'markerStyle', '.', 'markerColor', 'k');

oCfr = subplot(4,4,11);
imagesc(trialTimes, 1:size(spkMtx,2), corrCisMtx');
title('Odor C');
cdCorr = subplot(4,4,12);
corrScatPlot(corrCisMtx(:),corrDisMtx(:), 'C','D', 'markerStyle', '.', 'markerColor', 'k');

oDfr = subplot(4,4,16);
imagesc(trialTimes, 1:size(spkMtx,2), corrDisMtx');
title('Odor D');

frMapSPs = [oAfr, oBfr, oCfr, oDfr];
curCL = cell2mat(get(frMapSPs, 'clim'));
set(frMapSPs, 'clim', [min(curCL(:,1)), max(curCL(:,2))*.5], 'ydir', 'normal');
corrSPs = [abCorr, acCorr, adCorr, bcCorr, bdCorr, cdCorr];
curXL = cell2mat(get(corrSPs, 'xlim'));
curYL = cell2mat(get(corrSPs, 'ylim'));
set(corrSPs, 'xlim', [0 max(curXL(:,2))], 'ylim', [0 max(curYL(:,2))]);

axesHandles = findobj(get(gcf,'Children'), 'flat','Type','axes');
axis(axesHandles,'square');

annotation('textbox', 'position', [0.5 0.935 0.5 0.05], 'String', ['\bf\fontsize{10}' sprintf('Bin = %i ms; Step = %i ms', binSize, dsRate)],...
    'linestyle', 'none', 'horizontalalignment', 'right');
annotation('textbox', 'position', [0.025 0.935 0.5 0.05], 'String', '\bf\fontsize{14}Priors',...
    'linestyle', 'none', 'horizontalalignment', 'left');
curDir = cd;
annotation('textbox', 'position', [0.025 0.025 0.7 0.05], 'String', curDir,...
    'linestyle', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
colormap jet
orient(gcf, 'tall');
orient(gcf, 'landscape');
drawnow;


[~, aPriorAllPostRaw] = CalculatePostProb(corrAisMtx, nonAIStrials, binSize);
[~, bPriorAllPostRaw] = CalculatePostProb(corrBisMtx, nonAIStrials, binSize);
[~, cPriorAllPostRaw] = CalculatePostProb(corrCisMtx, nonAIStrials, binSize);
[~, dPriorAllPostRaw] = CalculatePostProb(corrDisMtx, nonAIStrials, binSize);
nonFISaLog = perfLog & inSeqLog & odorAlog & ~fullInSeqLog;
nonFISbLog = perfLog & inSeqLog & odorBlog & ~fullInSeqLog;
nonFIScLog = perfLog & inSeqLog & odorClog & ~fullInSeqLog;
nonFISdLog = perfLog & inSeqLog & odorDlog & ~fullInSeqLog;

allNonFISpost = [aPriorAllPostRaw, bPriorAllPostRaw, cPriorAllPostRaw, dPriorAllPostRaw];
for c = 1:size(allNonFISpost,1)
    allNonFISpost(c,:,:) = allNonFISpost(c,:,:)./sum(allNonFISpost(c,:,:));
end

aPriorAllPostNorm = allNonFISpost(:,1:size(aPriorAllPostRaw,1),:);
bPriorAllPostNorm = allNonFISpost(:,size(aPriorAllPostRaw,1)+1:size(aPriorAllPostRaw,1)*2,:);
cPriorAllPostNorm = allNonFISpost(:,size(aPriorAllPostRaw,1)*2+1:size(aPriorAllPostRaw,1)*3,:);
dPriorAllPostNorm = allNonFISpost(:,size(aPriorAllPostRaw,1)*3+1:end,:);

figure;
cAx = nan(4,4,2);
for prior = 1:4
    switch prior
        case 1
            curPrior = aPriorAllPostNorm;
        case 2
            curPrior = bPriorAllPostNorm;
        case 3
            curPrior = cPriorAllPostNorm;
        case 4
            curPrior = dPriorAllPostNorm;
    end
    for post = 1:4
        curPostLog = nonAISodors == post;
        curDecode = curPrior(:,:,curPostLog);
        subplot(4,4,sub2ind([4 4], post, prior))        
        cAx(prior,post,:) = PlotPostMtx(trialTimes, curDecode, sprintf('Prior%i; Decode%i', prior, post));
    end
end
cAxs = [min(min(cAx(:,:,1))), max(max(cAx(:,:,2)))*.5];

annotation('textbox', 'position', [0.5 0.935 0.5 0.05], 'String', ['\bf\fontsize{10}' sprintf('Bin = %i ms; Step = %i ms', binSize, dsRate)],...
    'linestyle', 'none', 'horizontalalignment', 'right');
annotation('textbox', 'position', [0.025 0.935 0.5 0.05], 'String', '\bf\fontsize{14}Decoding Trial and Odor Across Odors',...
    'linestyle', 'none', 'horizontalalignment', 'left');
curDir = cd;
annotation('textbox', 'position', [0.025 0.025 0.7 0.05], 'String', curDir,...
    'linestyle', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
colormap jet
axesHandles = findobj(get(gcf,'Children'), 'flat','Type','axes');
axis(axesHandles,'square')
set(axesHandles, 'clim', cAxs);
orient(gcf, 'tall');
orient(gcf, 'landscape');
drawnow;
            
%%
figure
corrISmtx = mean(spkMtx(:,:,perfLog & fullInSeqLog),3);                         % All InSeq Correct Trials
[post,~] = CalculatePostProb(corrISmtx, spkMtx(:,:,perfLog & inSeqLog & ~fullInSeqLog), binSize);
PlotPostMtx(trialTimes, post, 'InSeq Correct Trials');

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

%% Compare Trial Before with Trial After OutSeq Trial
% nonAIStrials = spkMtx(:,:,perfLog & inSeqLog & ~fullInSeqLog);
% nonAISodors = [nonAIStrialData.Odor];
% nonAIStrialNums = [nonAIStrialData.TrialNum];
% 
% preOStrial = false(size(nonAIStrialNums));
% postOStrial = false(size(nonAIStrialNums));
% for trl = 1:length(nonAIStrialNums)
%     if nonAIStrialNums(trl)~=1 && nonAIStrialNums(trl) ~= length(trialPeriodTD)
%         prevTrial = trialPeriodTD(nonAIStrialNums(trl)-1);
%         curTrl = trialPeriodTD(nonAIStrialNums(trl));
%         nextTrial = trialPeriodTD(nonAIStrialNums(trl)+1);
%         if curTrl.Position ~= 1
%             if curTrl.Position - prevTrial.Position == 1 && prevTrial.TranspositionDistance ~= 0
%                 postOStrial(trl) = true;
%                 preOStrial(trl-2) = true;
%             end
%         end
%    end
% end
% postPostOS = post(:,:,postOStrial);
% postPreOS = post(:,:,preOStrial);
% 
% figure;
% subplot(1,2,1)
% PlotPostMtx(trialTimes, postPreOS, 'Trial Before OutSeq');
% 
% subplot(1,2,2)
% PlotPostMtx(trialTimes, postPostOS, 'Trial After OutSeq');
% 
% annotation('textbox', 'position', [0.5 0.935 0.5 0.05], 'String', ['\bf\fontsize{10}' sprintf('Bin = %i ms; Step = %i ms', binSize, dsRate)],...
%     'linestyle', 'none', 'horizontalalignment', 'right');
% curDir = cd;
% annotation('textbox', 'position', [0.025 0.025 0.7 0.05], 'String', curDir,...
%     'linestyle', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
% colormap jet
% axesHandles = findobj(get(gcf,'Children'), 'flat','Type','axes');
% axis(axesHandles,'square')
% cLims = cell2mat(get(axesHandles, 'cLim'));
% set(axesHandles, 'cLim', [min(cLims(:,1)), max(cLims(:,2))]);
% orient(gcf, 'tall');
% orient(gcf, 'landscape');
% drawnow;
 
%% Comment in if you want to display each trial individually
% figure;
% for tr = 1:size(post,3)
%     imagesc(trialTimes, trialTimes, post(:,:,tr), [0 1]);
%     title(num2str(nonAISodors(tr)));
%     set(gca, 'ydir', 'normal');
%     pause(1);
% end
    
%%
function cAx = PlotPostMtx(trialTimes, postMtx, id)
imagesc(trialTimes, trialTimes, nanmean(postMtx,3));
title(id);
set(gca, 'ydir', 'normal')
hold on;
% line(trialTimes, trialTimes, 'linestyle', ':', 'color', 'w', 'linewidth', 2);
xlabel('True Time');
ylabel('Decoded Time');
drawnow
cAx = get(gca, 'clim');
end


%%
function [postNorm, postRaw] = CalculatePostProb(likely, obsv, binSize)
postNorm = nan(size(obsv,1), size(likely,1), size(obsv,3));
postRaw = nan(size(obsv,1), size(likely,1), size(obsv,3));
for trl = 1:size(obsv,3)
    for t = 1:size(obsv,1)
        p = nan(size(likely));
        curPopVect = obsv(t,:,trl)./(1000/binSize);
        curPopFact = factorial(curPopVect);
        for u = 1:size(likely,2)
            curAvgUniFR = likely(:,u);
            p(:,u) = (((binSize/1000).*curAvgUniFR).^curPopVect(u))./curPopFact(u);
        end        
        pp = prod(p,2);
        ee = exp(-(binSize/1000*sum(likely,2)));
        tempPost = pp.*ee;
        postRaw(t,:,trl) = tempPost;
        postNorm(t,:,trl) = tempPost./sum(tempPost);
    end
end
end