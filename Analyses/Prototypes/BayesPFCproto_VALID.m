%%
clc
clear all

%% Runtime variables
binSize = 50;
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

withdrawTime = 1.1;
%% Extract Behavioral Periods
% Taking 1/2 the binSize on either end to get rid of edge effects.
trialPeriodTD = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [-0.4-(binSize/2/1000) withdrawTime+0.5+(binSize/2/1000)], 'PokeIn');
% trialPeriodTD = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [-1.5-(binSize/2/1000) 0.5+(binSize/2/1000)], 'PokeOut');
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

%% Create Standard Likelihood Vector
fullyISlikelyCell = repmat({nan(size(spikeMatrix,1), size(spikeMatrix,2))}, 4,1);
odorCell = cell(4,1);
for p = 1:4
    fullyISlikelyCell{p} = mean(spikeMatrix(:,:,inSeqSeqs(p,:)),3);
    odorCell{p} = ones(size(fullyISlikelyCell{p},1),1)*p;
end
fullyISlikely = cell2mat(fullyISlikelyCell);
timeIndexes = repmat(trialTime,[4,1]);
odorIndexes = cell2mat(odorCell);

xTickNdxs = find(timeIndexes==0 | timeIndexes==min(timeIndexes) | timeIndexes==min(trialTime(trialTime>=withdrawTime)));
odorDivs = find(timeIndexes==min(timeIndexes));
odorTargs = find(timeIndexes==min(trialTime(trialTime>=withdrawTime)));
odorStarts = find(timeIndexes==0);

%% 
figure
subplot(4,5,1:4);
aPost = CalcStaticBayesPost(fullyISlikely, spikeMatrix(:,:,[trialPeriodTD.Odor]==1 & ~fullInSeqLog & perfLog & inSeqLog), binSize);
imagesc(1:length(timeIndexes), trialTime, nanmean(aPost,3),[0 0.0025]);
hold on;
set(gca, 'ydir', 'normal', 'xtick', xTickNdxs, 'xticklabel', repmat([-0.4 0 withdrawTime],[1,4]));
for p = 1:4
    plot([odorDivs(p) odorDivs(p)], get(gca, 'ylim'), '-w', 'linewidth', 2);
    plot([odorTargs(p) odorTargs(p)], get(gca, 'ylim'), ':w', 'linewidth', 1.5);
    plot([odorStarts(p) odorStarts(p)], get(gca, 'ylim'), '--w', 'linewidth', 1.5);
end
aDecode = DecodeBayesPost(aPost, odorIndexes);
aDecodePrcnts = nan(size(aDecode,1),4);
for t = 1:size(aDecode,1)
    aDecodePrcnts(t,1) = sum(aDecode(t,:)==1)/sum(~isnan(aDecode(t,:)));
    aDecodePrcnts(t,2) = sum(aDecode(t,:)==2)/sum(~isnan(aDecode(t,:)));
    aDecodePrcnts(t,3) = sum(aDecode(t,:)==3)/sum(~isnan(aDecode(t,:)));
    aDecodePrcnts(t,4) = sum(aDecode(t,:)==4)/sum(~isnan(aDecode(t,:)));
end
subplot(4,5,5);
curBar = barh(trialTime,aDecodePrcnts, 1, 'stacked');
curBar(1).FaceColor = [44/255 168/255 224/255];
curBar(2).FaceColor = [154/255 133/255 122/255]; 
curBar(3).FaceColor = [9/255 161/255 74/255];
curBar(4).FaceColor = [128/255 66/255 151/255];
set(gca, 'xlim', [0 1], 'xticklabel', []);
hold on
plot(get(gca, 'xlim'), [0 0], '--w');
plot(get(gca, 'xlim'), [withdrawTime withdrawTime], ':w');
box off
legend('A', 'B', 'C', 'D');
drawnow;
subplot(4,5,6:9);
bPost = CalcStaticBayesPost(fullyISlikely, spikeMatrix(:,:,[trialPeriodTD.Odor]==2 & ~fullInSeqLog & perfLog & inSeqLog), binSize);
imagesc(1:length(timeIndexes), trialTime, nanmean(bPost,3),[0 0.0025]);
hold on;
set(gca, 'ydir', 'normal', 'xtick', xTickNdxs, 'xticklabel', repmat([-0.4 0 withdrawTime],[1,4]));
for p = 1:4
    plot([odorDivs(p) odorDivs(p)], get(gca, 'ylim'), '-w', 'linewidth', 2);
    plot([odorTargs(p) odorTargs(p)], get(gca, 'ylim'), ':w', 'linewidth', 1.5);
    plot([odorStarts(p) odorStarts(p)], get(gca, 'ylim'), '--w', 'linewidth', 1.5);
end
bDecode = DecodeBayesPost(bPost, odorIndexes);
bDecodePrcnts = nan(size(bDecode,1),4);
for t = 1:size(aDecode,1)
    bDecodePrcnts(t,1) = sum(bDecode(t,:)==1)/sum(~isnan(bDecode(t,:)));
    bDecodePrcnts(t,2) = sum(bDecode(t,:)==2)/sum(~isnan(bDecode(t,:)));
    bDecodePrcnts(t,3) = sum(bDecode(t,:)==3)/sum(~isnan(bDecode(t,:)));
    bDecodePrcnts(t,4) = sum(bDecode(t,:)==4)/sum(~isnan(bDecode(t,:)));
end
subplot(4,5,10);
curBar = barh(trialTime,bDecodePrcnts, 1, 'stacked');
curBar(1).FaceColor = [44/255 168/255 224/255];
curBar(2).FaceColor = [154/255 133/255 122/255]; 
curBar(3).FaceColor = [9/255 161/255 74/255];
curBar(4).FaceColor = [128/255 66/255 151/255];
set(gca, 'xlim', [0 1], 'xticklabel', []);
hold on
plot(get(gca, 'xlim'), [0 0], '--w');
plot(get(gca, 'xlim'), [withdrawTime withdrawTime], ':w');
box off
drawnow;
subplot(4,5,11:14);
cPost = CalcStaticBayesPost(fullyISlikely, spikeMatrix(:,:,[trialPeriodTD.Odor]==3 & ~fullInSeqLog & perfLog & inSeqLog), binSize);
imagesc(1:length(timeIndexes), trialTime, nanmean(cPost,3),[0 0.0025]);
hold on;
set(gca, 'ydir', 'normal', 'xtick', xTickNdxs, 'xticklabel', repmat([-0.4 0 withdrawTime],[1,4]));
for p = 1:4
    plot([odorDivs(p) odorDivs(p)], get(gca, 'ylim'), '-w', 'linewidth', 2);
    plot([odorTargs(p) odorTargs(p)], get(gca, 'ylim'), ':w', 'linewidth', 1.5);
    plot([odorStarts(p) odorStarts(p)], get(gca, 'ylim'), '--w', 'linewidth', 1.5);
end
cDecode = DecodeBayesPost(cPost, odorIndexes);
cDecodePrcnts = nan(size(cDecode,1),4);
for t = 1:size(aDecode,1)
    cDecodePrcnts(t,1) = sum(cDecode(t,:)==1)/sum(~isnan(cDecode(t,:)));
    cDecodePrcnts(t,2) = sum(cDecode(t,:)==2)/sum(~isnan(cDecode(t,:)));
    cDecodePrcnts(t,3) = sum(cDecode(t,:)==3)/sum(~isnan(cDecode(t,:)));
    cDecodePrcnts(t,4) = sum(cDecode(t,:)==4)/sum(~isnan(cDecode(t,:)));
end
subplot(4,5,15);
curBar = barh(trialTime,cDecodePrcnts, 1, 'stacked');
curBar(1).FaceColor = [44/255 168/255 224/255];
curBar(2).FaceColor = [154/255 133/255 122/255]; 
curBar(3).FaceColor = [9/255 161/255 74/255];
curBar(4).FaceColor = [128/255 66/255 151/255];
set(gca, 'xlim', [0 1], 'xticklabel', []);
hold on
plot(get(gca, 'xlim'), [0 0], '--w');
plot(get(gca, 'xlim'), [withdrawTime withdrawTime], ':w');
box off
drawnow;
subplot(4,5,16:19);
dPost = CalcStaticBayesPost(fullyISlikely, spikeMatrix(:,:,[trialPeriodTD.Odor]==4 & ~fullInSeqLog & perfLog & inSeqLog), binSize);
imagesc(1:length(timeIndexes), trialTime, nanmean(dPost,3),[0 0.0025]);
hold on;
set(gca, 'ydir', 'normal', 'xtick', xTickNdxs, 'xticklabel', repmat([-0.4 0 withdrawTime],[1,4]));
for p = 1:4
    plot([odorDivs(p) odorDivs(p)], get(gca, 'ylim'), '-w', 'linewidth', 2);
    plot([odorTargs(p) odorTargs(p)], get(gca, 'ylim'), ':w', 'linewidth', 1.5);
    plot([odorStarts(p) odorStarts(p)], get(gca, 'ylim'), '--w', 'linewidth', 1.5);
end
dDecode = DecodeBayesPost(dPost, odorIndexes);
dDecodePrcnts = nan(size(dDecode,1),4);
for t = 1:size(aDecode,1)
    dDecodePrcnts(t,1) = sum(dDecode(t,:)==1)/sum(~isnan(dDecode(t,:)));
    dDecodePrcnts(t,2) = sum(dDecode(t,:)==2)/sum(~isnan(dDecode(t,:)));
    dDecodePrcnts(t,3) = sum(dDecode(t,:)==3)/sum(~isnan(dDecode(t,:)));
    dDecodePrcnts(t,4) = sum(dDecode(t,:)==4)/sum(~isnan(dDecode(t,:)));
end
subplot(4,5,20);
curBar = barh(trialTime,dDecodePrcnts, 1, 'stacked');
curBar(1).FaceColor = [44/255 168/255 224/255];
curBar(2).FaceColor = [154/255 133/255 122/255]; 
curBar(3).FaceColor = [9/255 161/255 74/255];
curBar(4).FaceColor = [128/255 66/255 151/255];
set(gca, 'xlim', [0 1], 'xticklabel', []);
hold on
plot(get(gca, 'xlim'), [0 0], '--w');
plot(get(gca, 'xlim'), [withdrawTime withdrawTime], ':w');
box off
drawnow;

colormap jet
%% 
function decode = DecodeBayesPost(post, id)
% Assumes post is in the structure of ObservTime X LikelyTime X Trial
decode = nan(size(post,1),size(post,3));
for o = 1:size(post,3)
    for d = 1:size(post,1)
        if ~isnan(post(d,1,o))
            decode(d,o) = id(find(post(d,:,o)==max(post(d,:,o)),1,'first'));
        end
    end
end
end
        
%%
function post = CalcStaticBayesPost(likely, obsv, binSize)
post = nan(size(obsv,1), size(likely,1), size(obsv,3));
for trl = 1:size(obsv,3)
    for t = 1:size(obsv,1)
        p = nan(size(likely));
        curPopVect = obsv(t,:,trl)*(binSize/1000);
        curPopFact = factorial(curPopVect);
        for u = 1:size(likely,2)
            curAvgUniFR = likely(:,u);
            p(:,u) = (((binSize/1000).*curAvgUniFR).^curPopVect(u))./curPopFact(u);
        end        
        pp = prod(p,2);
        ee = exp(-((binSize/1000)*sum(likely,2)));
        tempPost = pp.*ee;
        post(t,:,trl) = tempPost./sum(tempPost);
    end
end
end