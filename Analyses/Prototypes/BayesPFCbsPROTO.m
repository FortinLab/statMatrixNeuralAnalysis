%% Runtime variables
binSize = 200;
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
trialPeriodTD = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [-0.5 1.5], 'PokeIn');
% trialPeriodTD = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [-1.5 0.5], 'PokeOut');
trialEnsemble = ExtractTrialData_SM(trialPeriodTD, ensembleMatrix(:,2:end)); %#ok<*NODEF>
trialEnsembleMtx = cell2mat(reshape(trialEnsemble, [1 1 length(trialEnsemble)]));

trialTimes = behavMatrix(trialPeriodTD(1).TrialLogVect,1) - behavMatrix(trialPeriodTD(1).PokeInIndex,1);
% trialTimes = behavMatrix(trialPeriodTD(1).TrialLogVect,1) - behavMatrix(trialPeriodTD(1).PokeOutIndex,1);


%% Bin the spiking data
% First convolve the entire trialEnsembleMtx with a square to bin the
% spikes
binnedEnsembleMtx = nan(size(trialEnsembleMtx));
for t = 1:size(trialEnsembleMtx,3)
    for u = 1:size(trialEnsembleMtx,2)
        binnedEnsembleMtx(:,u,t) = conv(trialEnsembleMtx(:,u,t), ones(1,binSize)./(binSize/1000), 'same');
    end
end
% Now downsample the binned matrix
dsVect = downsample(1:size(binnedEnsembleMtx,1), dsRate);
spikeMatrix = binnedEnsembleMtx(dsVect,:,:);
trialTime = trialTimes(dsVect);

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

%% 
uniFRthreshLog = max(mean(spikeMatrix,3))<1;
spkMtx = spikeMatrix;
% spkMtx(:,uniFRthreshLog,:) = [];
% goodUniNames = {ensembleUnitSummaries(~uniFRthreshLog).UnitName};
%% Calculate Decoded Time Per Trial Type
figure;
corrAisMtx = mean(spkMtx(:,:,perfLog & fullInSeqLog & odorAlog),3);             % All A InSeq Correct Trials
[postA] = CalculatePostProb(corrAisMtx, spkMtx(:,:,perfLog & inSeqLog & odorAlog & ~fullInSeqLog), binSize);
subplot(2,2,1);
aCaxis = PlotPostMtx(trialTimes, postA, 'Odor A');

corrBisMtx = mean(spkMtx(:,:,perfLog & fullInSeqLog & odorBlog),3);             % All B InSeq Correct Trials
[postB] = CalculatePostProb(corrBisMtx, spkMtx(:,:,perfLog & inSeqLog & odorBlog & ~fullInSeqLog), binSize);
subplot(2,2,2);
bCaxis = PlotPostMtx(trialTimes, postB, 'Odor B');

corrCisMtx = mean(spkMtx(:,:,perfLog & fullInSeqLog & odorClog),3);             % All C InSeq Correct Trials
[postC] = CalculatePostProb(corrCisMtx, spkMtx(:,:,perfLog & inSeqLog & odorClog & ~fullInSeqLog), binSize);
subplot(2,2,3);
cCaxis = PlotPostMtx(trialTimes, postC, 'Odor C');

corrDisMtx = mean(spkMtx(:,:,perfLog & fullInSeqLog & odorDlog),3);             % All D InSeq Correct Trials
[postD] = CalculatePostProb(corrDisMtx, spkMtx(:,:,perfLog & inSeqLog & odorDlog & ~fullInSeqLog), binSize);
subplot(2,2,4);
dCaxis = PlotPostMtx(trialTimes, postD, 'Odor D');

cAx = [min([aCaxis, bCaxis, cCaxis, dCaxis]), max([aCaxis, bCaxis, cCaxis, dCaxis])];

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

%%
figure
corrISmtx = mean(spkMtx(:,:,perfLog & fullInSeqLog),3);                         % All InSeq Correct Trials
[post] = CalculatePostProb(corrISmtx, spkMtx(:,:,perfLog & inSeqLog & ~fullInSeqLog), binSize);
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
%%
function cAx = PlotPostMtx(trialTimes, postMtx, id)
imagesc(trialTimes, trialTimes, nanmean(postMtx,3)');
set(gca, 'ydir', 'normal')
title(id);
xlabel('True Time');
ylabel('Decoded Time');
drawnow
cAx = get(gca, 'clim');
end


%%
function [postMtx] = CalculatePostProb(meanFR, trialFR, binSize)
propVect = CalculateProportionalConstant(meanFR);
post = nan(size(trialFR,1), size(trialFR,1), size(trialFR,3));
% 
% for trl = 1:size(trialFR,3)
%     for i1 = 1:size(trialFR,1)
%         curPopVect = trialFR(i1,:,trl);
%         curPopVectFact = factorial(curPopVect);
%         for i2 = 1:size(trialFR,1)
%             curMeanFR = meanFR(i2,:);
%             condProbSpkPerUni = nan(size(curMeanFR));
%             for uni = 1:size(trialFR,2)
%                 condProbSpkPerUni(uni) = (((binSize/1000 * curMeanFR(uni))^curPopVect(uni))/curPopVectFact(uni)) * exp(-(binSize/1000*curMeanFR(uni)));
%             end
%             post(i1,i2,trl) = propVect(i1)*prod(condProbSpkPerUni);
%         end
%         post(i1,:,trl) = post(i1,:,trl)./max(post(i1,:,trl));
%     end
% end
% postMtx = post;

postPE = nan(size(trialFR,1), size(trialFR,1), size(trialFR,3));
for trl = 1:size(trialFR,3)
    for t = 1:size(meanFR,1)
        p = nan(size(meanFR));
%         e = nan(size(meanFR));
        curPopVect = trialFR(t,:,trl);
        curPopFact = factorial(curPopVect);
        for u = 1:size(meanFR,2)
            curAvgUniFR = meanFR(:,u);
            p(:,u) = (((binSize/1000)*curAvgUniFR).^curPopVect(u))./curPopFact(u);
%             e(:,u) = exp(-(binSize/1000*curAvgUniFR));
        end        
        pp = prod(p,2);
%         ee = sum(e,2);
        ee = exp(-(binSize/1000*sum(meanFR,2)));
        tempPost = propVect.*pp.*ee;
        postPE(t,:,trl) = tempPost./max(tempPost);
    end
end
postMtx = postPE;

end

%%
function [propConst] = CalculateProportionalConstant(rateMtx)
% propConst = nan(size(rateMtx));
% for u = 1:size(rateMtx,2)
%     propConstMtx(:,u) = 1./(rateMtx(:,u).*sum(rateMtx(:,u)~=0));
% end
% propConst(isinf(propConstMtx)) = 0;

% sum(rateMtx'.*(1./sum(rateMtx')))
propConst = 1./sum(rateMtx,2);
end
