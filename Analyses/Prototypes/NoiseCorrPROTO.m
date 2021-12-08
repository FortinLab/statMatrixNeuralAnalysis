function NoiseCorrPROTO
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
preTrialPeriodTD = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [-0.5 0], 'PokeIn');
preTrialEnsemble = ExtractTrialData_SM(preTrialPeriodTD, ensembleMatrix(:,2:end)); %#ok<*NODEF>
preTrialPopVect = cell2mat(cellfun(@(a)sum(a,1), preTrialEnsemble, 'uniformoutput', 0)');

rlyTrialPeriodTD = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [0 0.5], 'PokeIn');
rlyTrialEnsemble = ExtractTrialData_SM(rlyTrialPeriodTD, ensembleMatrix(:,2:end));
rlyTrialPopVect = cell2mat(cellfun(@(a)sum(a,1), rlyTrialEnsemble, 'uniformoutput', 0)');

latTrialPeriodTD = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [-0.5 0], 'PokeOut');
latTrialEnsemble = ExtractTrialData_SM(latTrialPeriodTD, ensembleMatrix(:,2:end));
latTrialPopVect = cell2mat(cellfun(@(a)sum(a,1), latTrialEnsemble, 'uniformoutput', 0)');

pstTrialPeriodTD = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [0 0.5], 'PokeOut');
pstTrialEnsemble = ExtractTrialData_SM(pstTrialPeriodTD, ensembleMatrix(:,2:end));
pstTrialPopVect = cell2mat(cellfun(@(a)sum(a,1), pstTrialEnsemble, 'uniformoutput', 0)');

allTrialPopVectZ = zscore([preTrialPopVect; rlyTrialPopVect; latTrialPopVect;pstTrialPopVect],0,1);

preTrialPopVectZall = allTrialPopVectZ(1:size(allTrialPopVectZ,1)/4,:);
rlyTrialPopVectZall = allTrialPopVectZ(size(allTrialPopVectZ,1)/4+1:size(allTrialPopVectZ,1)/2,:);
latTrialPopVectZall = allTrialPopVectZ(size(allTrialPopVectZ,1)/2+1:size(allTrialPopVectZ,1)*.75,:);
pstTrialPopVectZall = allTrialPopVectZ(size(allTrialPopVectZ,1)*.75+1:end,:);

%% Create Session Trial ID Logical Vectors
% Using preTrialPeriodTD here arbitrarilly, any of the TD variables above
% would work and they're all identical.
perfLog = [preTrialPeriodTD.Performance];
inSeqLog = [preTrialPeriodTD.TranspositionDistance]==0;
outSeqLog = [preTrialPeriodTD.TranspositionDistance]~=0 & abs([preTrialPeriodTD.TranspositionDistance])<10;
odorAlog = [preTrialPeriodTD.Odor] == 1;
odorBlog = [preTrialPeriodTD.Odor] == 2;
odorClog = [preTrialPeriodTD.Odor] == 3;
odorDlog = [preTrialPeriodTD.Odor] == 4;

fullInSeqSeqsStart = find(conv([preTrialPeriodTD.Odor], 1:4, 'valid')==20 & conv([preTrialPeriodTD.Position], 1:4, 'valid')==20 & conv([preTrialPeriodTD.Performance]*1, ones(1,4), 'valid')==4);
inSeqSeqs = nan(3,length(fullInSeqSeqsStart));
for iS = 1:length(fullInSeqSeqsStart)
    inSeqSeqs(1,iS) = fullInSeqSeqsStart(iS);
    inSeqSeqs(2,iS) = fullInSeqSeqsStart(iS) + 1;
    inSeqSeqs(3,iS) = fullInSeqSeqsStart(iS) + 2;
    inSeqSeqs(4,iS) = fullInSeqSeqsStart(iS) + 3;
end
fullInSeqLog = false(1,length(preTrialPeriodTD));
fullInSeqLog(inSeqSeqs(:)) = true;

%% Conv test/verification
% anss = nan(4);
% for a = 1:4
%     for b = 1:4
%         temp = 1:4;
%         temp(b) = a;
%         anss(a,b) = conv(temp, 1:4, 'valid');
%     end
% end
% anss
%%
% [corrTrlsPop, corrTrlsNoise] = ExaminePopAndNoiseCorrs(preTrialPopVectZall, rlyTrialPopVectZall, latTrialPopVectZall, pstTrialPopVectZall,...
%     perfLog & inSeqLog, 'ISC Trials');
% PlotPeriodDiffNoiseCorrs(corrTrlsNoise, 'ISC Trials');
% 
% [inCorrTrlsPop, inCorrTrlsNoise] = ExaminePopAndNoiseCorrs(preTrialPopVectZall, rlyTrialPopVectZall, latTrialPopVectZall, pstTrialPopVectZall,...
%     ~perfLog & inSeqLog, 'ISIc Trials');
% PlotPeriodDiffNoiseCorrs(inCorrTrlsNoise, 'ISIc Trials');
% 
% [aTrialsCorrPopISC, aTrialsCorrNoiseISC] = ExaminePopAndNoiseCorrs(preTrialPopVectZall, rlyTrialPopVectZall, latTrialPopVectZall, pstTrialPopVectZall,...
%     perfLog & odorAlog & inSeqLog, 'Odor A ISC');
% PlotPeriodDiffNoiseCorrs(aTrialsCorrNoiseISC, 'Odor A ISC');
% 
% [bTrialsCorrPopISC, bTrialsCorrNoiseISC] = ExaminePopAndNoiseCorrs(preTrialPopVectZall, rlyTrialPopVectZall, latTrialPopVectZall, pstTrialPopVectZall,...
%     perfLog & odorBlog & inSeqLog, 'Odor B ISC');
% PlotPeriodDiffNoiseCorrs(bTrialsCorrNoiseISC, 'Odor B ISC');
% 
% [cTrialsCorrPopISC, cTrialsCorrNoiseISC] = ExaminePopAndNoiseCorrs(preTrialPopVectZall, rlyTrialPopVectZall, latTrialPopVectZall, pstTrialPopVectZall,...
%     perfLog & odorClog & inSeqLog, 'Odor C ISC');
% PlotPeriodDiffNoiseCorrs(cTrialsCorrNoiseISC, 'Odor C ISC');
% 
% [dTrialsCorrPopISC, dTrialsCorrNoiseISC] = ExaminePopAndNoiseCorrs(preTrialPopVectZall, rlyTrialPopVectZall, latTrialPopVectZall, pstTrialPopVectZall,...
%     perfLog & odorDlog & inSeqLog, 'Odor D ISC');
% PlotPeriodDiffNoiseCorrs(dTrialsCorrNoiseISC, 'Odor D ISC');
% 
% [aTrialsCorrPopOSC, aTrialsCorrNoiseOSC] = ExaminePopAndNoiseCorrs(preTrialPopVectZall, rlyTrialPopVectZall, latTrialPopVectZall, pstTrialPopVectZall,...
%     perfLog & odorAlog & outSeqLog, 'Odor A OSC');
% PlotPeriodDiffNoiseCorrs(aTrialsCorrNoiseOSC, 'Odor A OSC');
% 
% [bTrialsCorrPopOSC, bTrialsCorrNoiseOSC] = ExaminePopAndNoiseCorrs(preTrialPopVectZall, rlyTrialPopVectZall, latTrialPopVectZall, pstTrialPopVectZall,...
%     perfLog & odorBlog & outSeqLog, 'Odor B OSC');
% PlotPeriodDiffNoiseCorrs(bTrialsCorrNoiseOSC, 'Odor B OSC');
% 
% [cTrialsCorrPopOSC, cTrialsCorrNoiseOSC] = ExaminePopAndNoiseCorrs(preTrialPopVectZall, rlyTrialPopVectZall, latTrialPopVectZall, pstTrialPopVectZall,...
%     perfLog & odorClog & outSeqLog, 'Odor C OSC');
% PlotPeriodDiffNoiseCorrs(cTrialsCorrNoiseOSC, 'Odor C OSC');
% 
% [dTrialsCorrPopOSC, dTrialsCorrNoiseOSC] = ExaminePopAndNoiseCorrs(preTrialPopVectZall, rlyTrialPopVectZall, latTrialPopVectZall, pstTrialPopVectZall,...
%     perfLog & odorDlog & outSeqLog, 'Odor D OSC');
% PlotPeriodDiffNoiseCorrs(dTrialsCorrNoiseOSC, 'Odor D OSC');

%% Only Trials from Complete InSeq Seqs
% [aTrialsCorrPopCompISC, aTrialsCorrNoiseCompISC] = ExaminePopAndNoiseCorrs(preTrialPopVectZall, rlyTrialPopVectZall, latTrialPopVectZall, pstTrialPopVectZall,...
%     perfLog & odorAlog & fullInSeqLog, 'Odor A Complete Seqs');
% PlotPeriodDiffNoiseCorrs(aTrialsCorrNoiseCompISC, 'Odor A Complete Seqs');
% 
% [bTrialsCorrPopCompISC, bTrialsCorrNoiseCompISC] = ExaminePopAndNoiseCorrs(preTrialPopVectZall, rlyTrialPopVectZall, latTrialPopVectZall, pstTrialPopVectZall,...
%     perfLog & odorBlog & fullInSeqLog, 'Odor B Complete Seqs');
% PlotPeriodDiffNoiseCorrs(bTrialsCorrNoiseCompISC, 'Odor B Complete Seqs');
% 
% [cTrialsCorrPopCompISC, cTrialsCorrNoiseCompISC] = ExaminePopAndNoiseCorrs(preTrialPopVectZall, rlyTrialPopVectZall, latTrialPopVectZall, pstTrialPopVectZall,...
%     perfLog & odorClog & fullInSeqLog, 'Odor C Complete Seqs');
% PlotPeriodDiffNoiseCorrs(cTrialsCorrNoiseCompISC, 'Odor C Complete Seqs');
% 
% [dTrialsCorrPopCompISC, dTrialsCorrNoiseCompISC] = ExaminePopAndNoiseCorrs(preTrialPopVectZall, rlyTrialPopVectZall, latTrialPopVectZall, pstTrialPopVectZall,...
%     perfLog & odorDlog & fullInSeqLog, 'Odor D Complete Seqs');
% PlotPeriodDiffNoiseCorrs(dTrialsCorrNoiseCompISC, 'Odor D Complete Seqs');
% 
% PlotOdorDiffNoiseCorrs(aTrialsCorrNoiseCompISC, bTrialsCorrNoiseCompISC, cTrialsCorrNoiseCompISC, dTrialsCorrNoiseCompISC, 'Complete IS')

%% Run Permutations to derive 'significant' noise correlation values
% [~, aTrialsCorrNoiseCompISCz, ~, ~] = ExaminePopAndNoiseCorrsPerm(preTrialPopVectZall, rlyTrialPopVectZall, latTrialPopVectZall, pstTrialPopVectZall,...
%     perfLog & odorAlog & fullInSeqLog, 1000, 'Odor A Complete Seqs');
% 
% [~, bTrialsCorrNoiseCompISCz, ~, ~] = ExaminePopAndNoiseCorrsPerm(preTrialPopVectZall, rlyTrialPopVectZall, latTrialPopVectZall, pstTrialPopVectZall,...
%     perfLog & odorBlog & fullInSeqLog, 1000, 'Odor B Complete Seqs');
% 
% [~, cTrialsCorrNoiseCompISCz, ~, ~] = ExaminePopAndNoiseCorrsPerm(preTrialPopVectZall, rlyTrialPopVectZall, latTrialPopVectZall, pstTrialPopVectZall,...
%     perfLog & odorClog & fullInSeqLog, 1000, 'Odor C Complete Seqs');
% 
% [~, dTrialsCorrNoiseCompISCz, ~, ~] = ExaminePopAndNoiseCorrsPerm(preTrialPopVectZall, rlyTrialPopVectZall, latTrialPopVectZall, pstTrialPopVectZall,...
%     perfLog & odorDlog & fullInSeqLog, 1000, 'Odor D Complete Seqs');

[~, aTrialsCorrNoiseCompISCz, ~, ~] = ExaminePopAndNoiseCorrsPerm(preTrialPopVectZall, rlyTrialPopVectZall, latTrialPopVectZall, pstTrialPopVectZall,...
    perfLog & odorAlog & inSeqLog, 1000, 'Odor A ISC');

[~, bTrialsCorrNoiseCompISCz, ~, ~] = ExaminePopAndNoiseCorrsPerm(preTrialPopVectZall, rlyTrialPopVectZall, latTrialPopVectZall, pstTrialPopVectZall,...
    perfLog & odorBlog & inSeqLog, 1000, 'Odor B ISC');

[~, cTrialsCorrNoiseCompISCz, ~, ~] = ExaminePopAndNoiseCorrsPerm(preTrialPopVectZall, rlyTrialPopVectZall, latTrialPopVectZall, pstTrialPopVectZall,...
    perfLog & odorClog & inSeqLog, 1000, 'Odor C ISC');

[~, dTrialsCorrNoiseCompISCz, ~, ~] = ExaminePopAndNoiseCorrsPerm(preTrialPopVectZall, rlyTrialPopVectZall, latTrialPopVectZall, pstTrialPopVectZall,...
    perfLog & odorDlog & inSeqLog, 1000, 'Odor D ISC');

% PlotOdorDiffNoiseCorrs(aTrialsCorrNoiseCompISCz, bTrialsCorrNoiseCompISCz, cTrialsCorrNoiseCompISCz, dTrialsCorrNoiseCompISCz, 'Complete IS')


%%
% PlotOdorDiffNoiseCorrs(aTrialsCorrNoiseISC, bTrialsCorrNoiseISC, cTrialsCorrNoiseISC, dTrialsCorrNoiseISC, 'ISC')
% PlotOdorDiffNoiseCorrs(aTrialsCorrNoiseOSC, bTrialsCorrNoiseOSC, cTrialsCorrNoiseOSC, dTrialsCorrNoiseOSC, 'OSC')


end

function [popCorrsZ, noiseCorrsZ, popCorrsChance, noiseCorrsChance] = ExaminePopAndNoiseCorrsPerm(preTrial, rlyTrial, latTrial, pstTrial, sortVect, numPerms, sortVectID)
%% Chance up that random stream
randperm(floor(sum(clock)));
threshVal = 0;
%%
preSorted = preTrial(sortVect,:);
rlySorted = rlyTrial(sortVect,:);
latSorted = latTrial(sortVect,:);
pstSorted = pstTrial(sortVect,:);

noiseCorrs{1} = corr(preSorted);
popCorrs{1} = corr(preSorted');

noiseCorrs{2} = corr(rlySorted);
popCorrs{2} = corr(rlySorted');

noiseCorrs{3} = corr(latSorted);
popCorrs{3} = corr(latSorted');

noiseCorrs{4} = corr(pstSorted);
popCorrs{4} = corr(pstSorted');

noiseCorrsZ = repmat({nan(size(noiseCorrs{1}))}, 1, 4);
popCorrsZ = repmat({nan(size(popCorrs{1}))}, 1, 4);
noiseCorrsChance = repmat({nan(size(noiseCorrs{1},1), size(noiseCorrs{1},2), numPerms)},1,4);
popCorrsChance = repmat({nan(size(popCorrs{1},1), size(popCorrs{1},2), numPerms)},1,4);
for perm = 1:numPerms
    tempPre = preSorted;
    tempRly = rlySorted;
    tempLat = latSorted;
    tempPst = pstSorted;
    
    preIDmtx = reshape(1:length(tempPre(:)), size(tempPre,1), size(tempPre,2));
    for c = 1:size(preIDmtx,2)
        preIDmtx(:,c) = preIDmtx(randperm(size(preIDmtx,1)),c);
    end
    noiseCorrsChance{1}(:,:,perm) = corr(tempPre(preIDmtx));
    popCorrsChance{1}(:,:,perm) = corr(tempPre(preIDmtx)');
    
    rlyIDmtx = reshape(1:length(tempRly(:)), size(tempRly,1), size(tempRly,2));
    for c = 1:size(rlyIDmtx,2)
        rlyIDmtx(:,c) = rlyIDmtx(randperm(size(rlyIDmtx,1)),c);
    end
    noiseCorrsChance{2}(:,:,perm) = corr(tempRly(rlyIDmtx));
    popCorrsChance{2}(:,:,perm) = corr(tempRly(rlyIDmtx)');
    
    latIDmtx = reshape(1:length(tempLat(:)), size(tempLat,1), size(tempLat,2));
    for c = 1:size(latIDmtx,2)
        latIDmtx(:,c) = latIDmtx(randperm(size(latIDmtx,1)),c);
    end
    noiseCorrsChance{3}(:,:,perm) = corr(tempLat(latIDmtx));
    popCorrsChance{3}(:,:,perm) = corr(tempLat(latIDmtx)');
    
    pstIDmtx = reshape(1:length(tempPst(:)), size(tempPst,1), size(tempPst,2));
    for c = 1:size(pstIDmtx,2)
        pstIDmtx(:,c) = pstIDmtx(randperm(size(pstIDmtx,1)),c);
    end
    noiseCorrsChance{4}(:,:,perm) = corr(tempPst(pstIDmtx));
    popCorrsChance{4}(:,:,perm) = corr(tempPst(pstIDmtx)');
end
for nIr = 1:size(noiseCorrsZ{1},1)
    for nIc = 1:size(noiseCorrsZ{1},2)
        tempZpre = zscore([reshape(noiseCorrsChance{1}(nIr,nIc,:), 1, numPerms), noiseCorrs{1}(nIr,nIc)]);
        noiseCorrsZ{1}(nIr,nIc) = tempZpre(end);
        
        tempZrly = zscore([reshape(noiseCorrsChance{2}(nIr,nIc,:), 1, numPerms), noiseCorrs{2}(nIr,nIc)]);
        noiseCorrsZ{2}(nIr,nIc) = tempZrly(end);
        
        tempZlat = zscore([reshape(noiseCorrsChance{3}(nIr,nIc,:), 1, numPerms), noiseCorrs{3}(nIr,nIc)]);
        noiseCorrsZ{3}(nIr,nIc) = tempZlat(end);
        
        tempZpst = zscore([reshape(noiseCorrsChance{4}(nIr,nIc,:), 1, numPerms), noiseCorrs{4}(nIr,nIc)]);
        noiseCorrsZ{4}(nIr,nIc) = tempZpst(end);
    end
end

for pIr = 1:size(popCorrsZ{1},1)
    for pIc = 1:size(popCorrsZ{1},2)
        tempZpre = zscore([reshape(popCorrsChance{1}(pIr,pIc,:), 1, numPerms), popCorrs{1}(pIr,pIc)]);
        popCorrsZ{1}(pIr,pIc) = tempZpre(end);
        
        tempZrly = zscore([reshape(popCorrsChance{2}(pIr,pIc,:), 1, numPerms), popCorrs{2}(pIr,pIc)]);
        popCorrsZ{2}(pIr,pIc) = tempZrly(end);
        
        tempZlat = zscore([reshape(popCorrsChance{3}(pIr,pIc,:), 1, numPerms), popCorrs{3}(pIr,pIc)]);
        popCorrsZ{3}(pIr,pIc) = tempZlat(end);
        
        tempZpst = zscore([reshape(popCorrsChance{4}(pIr,pIc,:), 1, numPerms), popCorrs{4}(pIr,pIc)]);
        popCorrsZ{4}(pIr,pIc) = tempZpst(end);
    end
end

maxAbsZNoise = max(max(abs(cell2mat(noiseCorrsZ))));
maxAbsZPop = max(max(abs(cell2mat(popCorrsZ))));

figure; 
subplot(4,4,1)
PlotCombinedMatrix(noiseCorrs{1}, noiseCorrsZ{1}, maxAbsZNoise);
title('Pre Trial')
xlabel('Z-Normed')
ylabel('Observed')

subplot(4,4,2)
PlotCombinedMatrix(noiseCorrsZ{2}, noiseCorrsZ{1}, maxAbsZNoise);
title('Pre vs Early')
xlabel('Pre');
ylabel('Early');

subplot(4,4,3)
PlotCombinedMatrix(noiseCorrsZ{3}, noiseCorrsZ{1}, maxAbsZNoise);
title('Pre vs Late')
xlabel('Pre');
ylabel('Late');

subplot(4,4,4)
PlotCombinedMatrix(noiseCorrsZ{4}, noiseCorrsZ{1}, maxAbsZNoise);
title('Pre vs Post')
xlabel('Pre');
ylabel('Post');

subplot(4,4,5)
PlotNoiseCorrChanceDist(noiseCorrsChance{1}, noiseCorrsChance{2}, noiseCorrsZ{1}, noiseCorrsZ{2});
% upperPre = ExtractUpperTriangle(noiseCorrsZ{1});
% upperEarly = ExtractUpperTriangle(noiseCorrsZ{2});
% log = (~isnan(upperPre) & ~isnan(upperEarly)) & (abs(upperPre)>=threshVal | abs(upperEarly)>=threshVal);
% corrScatPlot(upperPre(log),upperEarly(log), 'Pre','Early',[]);
title('Pre vs Early');

subplot(4,4,6)
PlotCombinedMatrix(noiseCorrs{2}, noiseCorrsZ{2}, maxAbsZNoise);
title('Early Trial')
xlabel('Z-Normed')
ylabel('Observed')

subplot(4,4,7)
PlotCombinedMatrix(noiseCorrsZ{3}, noiseCorrsZ{2}, maxAbsZNoise);
title('Early vs Late')
xlabel('Early')
ylabel('Late')

subplot(4,4,8)
PlotCombinedMatrix(noiseCorrsZ{4}, noiseCorrsZ{2}, maxAbsZNoise);
title('Early vs Post')
xlabel('Early')
ylabel('Post')

subplot(4,4,9)
PlotNoiseCorrChanceDist(noiseCorrsChance{1}, noiseCorrsChance{3}, noiseCorrsZ{1}, noiseCorrsZ{3});
% upperPre = ExtractUpperTriangle(noiseCorrsZ{1});
% upperLate = ExtractUpperTriangle(noiseCorrsZ{3});
% log = (~isnan(upperPre) & ~isnan(upperLate)) & (abs(upperPre)>=threshVal | abs(upperPre)>=threshVal);
% corrScatPlot(upperPre(log),upperLate(log), 'Pre','Late',[]);
title('Pre vs Late')

subplot(4,4,10)
PlotNoiseCorrChanceDist(noiseCorrsChance{2}, noiseCorrsChance{3}, noiseCorrsZ{2}, noiseCorrsZ{3});
% upperEarly = ExtractUpperTriangle(noiseCorrsZ{2});
% upperLate = ExtractUpperTriangle(noiseCorrsZ{3});
% log = (~isnan(upperEarly) & ~isnan(upperLate)) & (abs(upperEarly)>=threshVal | abs(upperLate)>=threshVal);
% corrScatPlot(upperEarly(log),upperLate(log), 'Early','Late',[]);
title('Early vs Late')

subplot(4,4,11)
PlotCombinedMatrix(noiseCorrs{3}, noiseCorrsZ{3}, maxAbsZNoise);
title('Late Trial')
xlabel('Z-Normed')
ylabel('Observed')

subplot(4,4,12)
PlotCombinedMatrix(noiseCorrsZ{4}, noiseCorrsZ{3}, maxAbsZNoise);
title('Late vs Post')
xlabel('Late')
ylabel('Post')

subplot(4,4,13)
PlotNoiseCorrChanceDist(noiseCorrsChance{1}, noiseCorrsChance{4}, noiseCorrsZ{1}, noiseCorrsZ{4});
% upperPre = ExtractUpperTriangle(noiseCorrsZ{1});
% upperPost = ExtractUpperTriangle(noiseCorrsZ{4});
% log = (~isnan(upperPre) & ~isnan(upperPost)) & (abs(upperPre)>=threshVal | abs(upperPost)>=threshVal);
% corrScatPlot(upperPre(log),upperPost(log), 'Pre','Post',[]);
title('Pre vs Post')

subplot(4,4,14)
PlotNoiseCorrChanceDist(noiseCorrsChance{2}, noiseCorrsChance{4}, noiseCorrsZ{2}, noiseCorrsZ{4});
% upperEarly = ExtractUpperTriangle(noiseCorrsZ{2});
% upperPost = ExtractUpperTriangle(noiseCorrsZ{4});
% log = (~isnan(upperEarly) & ~isnan(upperPost)) & (abs(upperEarly)>=threshVal | abs(upperPost)>=threshVal);
% corrScatPlot(upperEarly(log),upperPost(log), 'Early','Post',[]);
title('Early vs Post')

subplot(4,4,15)
PlotNoiseCorrChanceDist(noiseCorrsChance{3}, noiseCorrsChance{4}, noiseCorrsZ{3}, noiseCorrsZ{4});
% upperLate = ExtractUpperTriangle(noiseCorrsZ{3});
% upperPost = ExtractUpperTriangle(noiseCorrsZ{4});
% log = (~isnan(upperLate) & ~isnan(upperPost)) & (abs(upperLate)>=threshVal | abs(upperPost)>=threshVal);
% corrScatPlot(upperLate(log),upperPost(log), 'Late','Post',[]);
title('Late vs Post')

subplot(4,4,16)
PlotCombinedMatrix(noiseCorrs{4}, noiseCorrsZ{4}, maxAbsZNoise);
title('Post Trial')
xlabel('Z-Normed')
ylabel('Observed')
    
annotation('textbox', 'position', [0.025 0.935 0.7 0.05], 'String', ['\bf\fontsize{14}' sortVectID],...
    'linestyle', 'none', 'horizontalalignment', 'left');
curDir = cd;
annotation('textbox', 'position', [0.025 0.025 0.7 0.05], 'String', curDir,...
    'linestyle', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');

colormap jet
axesHandles = findobj(get(gcf,'Children'), 'flat','Type','axes');
axis(axesHandles,'square')
orient(gcf, 'tall');
orient(gcf, 'landscape');

% 
% figure; 
% subplot(4,4,3)
% PlotCombinedMatrix(popCorrs{1}, popCorrsZ{1}, maxAbsZPop);
% subplot(4,4,7)
% PlotCombinedMatrix(popCorrs{2}, popCorrsZ{2}, maxAbsZPop);
% subplot(4,4,11)
% PlotCombinedMatrix(popCorrs{3}, popCorrsZ{3}, maxAbsZPop);
% subplot(4,4,15)
% PlotCombinedMatrix(popCorrs{4}, popCorrsZ{4}, maxAbsZPop);
% 
% colormap jet
% 
% axesHandles = findobj(get(gcf,'Children'), 'flat','Type','axes');
% axis(axesHandles,'square')
% orient(gcf, 'tall');
% orient(gcf, 'landscape');

end

%%
function PlotNoiseCorrChanceDist(chanceX, chanceY, realX, realY)
[upperRealX, ~] = ExtractUpperTriangle(realX);
[upperRealY, ~] = ExtractUpperTriangle(realY);
log = ~isnan(upperRealX) & ~isnan(upperRealY);
realCorr = corr(upperRealX(log), upperRealY(log));
chanceCorr = nan(1,size(chanceX,3));
for perm = 1:size(chanceX,3)
    [upperX, ~] = ExtractUpperTriangle(chanceX(:,:,perm));
    [upperY, ~] = ExtractUpperTriangle(chanceY(:,:,perm));
    log = ~isnan(upperX) & ~isnan(upperY);
    chanceCorr(perm) = corr(upperX(log), upperY(log));
end
histogram(chanceCorr, -0.2:0.01:.2)
hold on;
line([realCorr, realCorr], get(gca, 'ylim'), 'color', 'r', 'linewidth', 2);
box off;
grid on;
set(gca, 'color', 'none');
temp = zscore([chanceCorr, realCorr]);
xlabel(sprintf('Z = %.02f', temp(end)))
end
%%
function PlotCombinedMatrix(realVals, zVals, cMax)
[~, realMtxNanD] = ExtractUpperTriangle(realVals);
[~, zMtxNanD] = ExtractUpperTriangle(zVals);
imagesc(nansum(reshape([fliplr(realMtxNanD), fliplr(zMtxNanD')], size(realMtxNanD,1), size(realMtxNanD,2),2),3), [-cMax cMax]);
end


%%
function PlotPeriodDiffNoiseCorrs(mtx, id)
figure;
subplot(4,4,1)
imagesc(mtx{1}, [-1 1]);
title('Pre Trial');

subplot(4,4,2)
preUpper = ExtractUpperTriangle(mtx{1});
rlyUpper = ExtractUpperTriangle(mtx{2});
corrScatPlot(preUpper(~isnan(preUpper) & ~isnan(rlyUpper)),rlyUpper(~isnan(preUpper) & ~isnan(rlyUpper)), 'Pre','Early',[]);

subplot(4,4,3)
preUpper = ExtractUpperTriangle(mtx{1});
latUpper = ExtractUpperTriangle(mtx{3});
corrScatPlot(preUpper(~isnan(preUpper) & ~isnan(latUpper)),latUpper(~isnan(preUpper) & ~isnan(latUpper)), 'Pre','Late',[]);

subplot(4,4,4)
preUpper = ExtractUpperTriangle(mtx{1});
pstUpper = ExtractUpperTriangle(mtx{4});
corrScatPlot(preUpper(~isnan(preUpper) & ~isnan(pstUpper)),pstUpper(~isnan(preUpper) & ~isnan(pstUpper)), 'Pre','Post',[]);

subplot(4,4,5)
imagesc(mtx{1} - mtx{2}, [-1 1]);
title('Pre - Early');

subplot(4,4,6)
imagesc(mtx{2}, [-1 1]);
title('Early Trial')

subplot(4,4,7)
rlyUpper = ExtractUpperTriangle(mtx{2});
latUpper = ExtractUpperTriangle(mtx{3});
corrScatPlot(rlyUpper(~isnan(rlyUpper) & ~isnan(latUpper)),latUpper(~isnan(rlyUpper) & ~isnan(latUpper)), 'Early','Late',[]);

subplot(4,4,8)
rlyUpper = ExtractUpperTriangle(mtx{2});
pstUpper = ExtractUpperTriangle(mtx{4});
corrScatPlot(rlyUpper(~isnan(rlyUpper) & ~isnan(pstUpper)),pstUpper(~isnan(rlyUpper) & ~isnan(pstUpper)), 'Early','Post',[]);

subplot(4,4,9)
imagesc(mtx{1}-mtx{3}, [-1 1]);
title('Pre - Late');

subplot(4,4,10)
imagesc(mtx{2}-mtx{3}, [-1 1]);
title('Early - Late');

subplot(4,4,11)
imagesc(mtx{3}, [-1 1]);
title('Late Trial');

subplot(4,4,12)
latUpper = ExtractUpperTriangle(mtx{3});
pstUpper = ExtractUpperTriangle(mtx{4});
corrScatPlot(latUpper(~isnan(latUpper) & ~isnan(pstUpper)),pstUpper(~isnan(latUpper) & ~isnan(pstUpper)), 'Early','Post',[]);

subplot(4,4,13)
imagesc(mtx{1}-mtx{4}, [-1 1]);
title('Pre - Post');

subplot(4,4,14)
imagesc(mtx{2}-mtx{4}, [-1 1]);
title('Early - Post');

subplot(4,4,15)
imagesc(mtx{3}-mtx{4}, [-1 1]);
title('Late - Post');

subplot(4,4,16)
imagesc(mtx{4}, [-1 1]);
title('Post Trial');

annotation('textbox', 'position', [0.025 0.935 0.7 0.05], 'String', ['\bf\fontsize{14}' id],...
    'linestyle', 'none', 'horizontalalignment', 'left');

curDir = cd;
annotation('textbox', 'position', [0.025 0.025 0.7 0.05], 'String', curDir,...
    'linestyle', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');

colormap jet

axesHandles = findobj(get(gcf,'Children'), 'flat','Type','axes');
axis(axesHandles,'square')
orient(gcf, 'tall');
orient(gcf, 'landscape');
end

%% 
function [upperMtxVect, lowerMtxNand] = ExtractUpperTriangle(mtx)
upperLog = triu(true(size(mtx,1)),1);
upperMtxVect = mtx(upperLog);
lowerMtxNand = mtx;
lowerMtxNand(~upperLog) = nan;
end

%% 
function PlotOdorDiffNoiseCorrs(aTrialsCorrNoise, bTrialsCorrNoise, cTrialsCorrNoise, dTrialsCorrNoise, id)
trlPrdIDs = [{'Pre-Trial'}, {'Early-Trial'}, {'Late-Trial'}, {'Post-Trial'}];
for trlPrd = 1:4
    currA = aTrialsCorrNoise{trlPrd};
    currB = bTrialsCorrNoise{trlPrd};
    currC = cTrialsCorrNoise{trlPrd};
    currD = dTrialsCorrNoise{trlPrd};    
    figure;
    abDiff = currA - currB;
    subplot(4,4,2)
    imagesc(abDiff, [-1 1]);
    title('A-B')
    subplot(4,4,5)
    upperA = ExtractUpperTriangle(currA);
    upperB = ExtractUpperTriangle(currB);
    corrScatPlot(upperA(~isnan(upperA) & ~isnan(upperB)),upperB(~isnan(upperA) & ~isnan(upperB)), 'A','B',[]);
    set(gca, 'xlim', [-1 1], 'ylim', [-1 1]);
%     histogram(abDiff(logical(triu(ones(size(abDiff)),1))), -1:0.01:1, 'edgecolor', 'none', 'facecolor', 'k');
    set(gca, 'color', 'none', 'box', 'off');
    grid on
    hold on
%     line([0 0], get(gca,'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
    title('A vs B')
    
    acDiff = currA - currC;
    subplot(4,4,3)
    imagesc(acDiff, [-1 1]);
    title('A-C')
    subplot(4,4,9)
    upperA = ExtractUpperTriangle(currA);
    upperC = ExtractUpperTriangle(currC);
    corrScatPlot(upperA(~isnan(upperA) & ~isnan(upperC)),upperC(~isnan(upperA) & ~isnan(upperC)), 'A','C',[]);
    set(gca, 'xlim', [-1 1], 'ylim', [-1 1]);
%     histogram(acDiff(logical(triu(ones(size(acDiff)),1))), -1:0.01:1, 'edgecolor', 'none', 'facecolor', 'k');
    set(gca, 'color', 'none', 'box', 'off');
    grid on
    hold on
%     line([0 0], get(gca,'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
    title('A vs C')
    
    adDiff = currA - currD;
    subplot(4,4,4)
    imagesc(adDiff, [-1 1]);
    title('A-D')
    subplot(4,4,13)
    upperA = ExtractUpperTriangle(currA);
    upperD = ExtractUpperTriangle(currD);
    corrScatPlot(upperA(~isnan(upperA) & ~isnan(upperD)),upperD(~isnan(upperA) & ~isnan(upperD)), 'A','D',[]);
    set(gca, 'xlim', [-1 1], 'ylim', [-1 1]);
%     histogram(adDiff(logical(triu(ones(size(adDiff)),1))), -1:0.01:1, 'edgecolor', 'none', 'facecolor', 'k');
    set(gca, 'color', 'none', 'box', 'off');
    grid on
    hold on
%     line([0 0], get(gca,'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
    title('A vs D')
    
    bcDiff = currB - currC;
    subplot(4,4,7)
    imagesc(bcDiff, [-1 1]);
    title('B-C')
    subplot(4,4,10)
    upperB = ExtractUpperTriangle(currB);
    upperC = ExtractUpperTriangle(currC);
    corrScatPlot(upperB(~isnan(upperB) & ~isnan(upperC)),upperC(~isnan(upperB) & ~isnan(upperC)), 'B','C',[]);
    set(gca, 'xlim', [-1 1], 'ylim', [-1 1]);
%     histogram(bcDiff(logical(triu(ones(size(bcDiff)),1))), -1:0.01:1, 'edgecolor', 'none', 'facecolor', 'k');
    set(gca, 'color', 'none', 'box', 'off');
    grid on
    hold on
%     line([0 0], get(gca,'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
    title('B vs C')
    
    bdDiff = currB - currD;
    subplot(4,4,8)
    imagesc(bdDiff, [-1 1]);
    title('B-D')
    subplot(4,4,14)
    upperB = ExtractUpperTriangle(currB);
    upperD = ExtractUpperTriangle(currD);
    corrScatPlot(upperB(~isnan(upperB) & ~isnan(upperD)),upperD(~isnan(upperB) & ~isnan(upperD)), 'B','D',[]);
    set(gca, 'xlim', [-1 1], 'ylim', [-1 1]);
%     histogram(bdDiff(logical(triu(ones(size(bdDiff)),1))), -1:0.01:1, 'edgecolor', 'none', 'facecolor', 'k');
    set(gca, 'color', 'none', 'box', 'off');
    grid on
    hold on
%     line([0 0], get(gca,'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
    title('B vs D')
    
    cdDiff = currC - currD;
    subplot(4,4,12)
    imagesc(cdDiff, [-1 1]);
    title('C-D')
    subplot(4,4,15)
    upperC = ExtractUpperTriangle(currC);
    upperD = ExtractUpperTriangle(currD);
    corrScatPlot(upperC(~isnan(upperC) & ~isnan(upperD)),upperD(~isnan(upperC) & ~isnan(upperD)), 'C','D',[]);
    set(gca, 'xlim', [-1 1], 'ylim', [-1 1]);
%     histogram(cdDiff(logical(triu(ones(size(cdDiff)),1))), -1:0.01:1, 'edgecolor', 'none', 'facecolor', 'k');
    set(gca, 'color', 'none', 'box', 'off');
    grid on
    hold on
%     line([0 0], get(gca,'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
    title('C vs D')
    
    annotation('textbox', 'position', [0.025 0.935 0.7 0.05], 'String', ['\bf\fontsize{14}' sprintf('%s %s', trlPrdIDs{trlPrd}, id)],...
        'linestyle', 'none', 'horizontalalignment', 'left');
    
    curDir = cd;
    annotation('textbox', 'position', [0.025 0.025 0.7 0.05], 'String', curDir,...
        'linestyle', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
    
    colormap jet

    axesHandles = findobj(get(gcf,'Children'), 'flat','Type','axes');
    axis(axesHandles,'square')
    orient(gcf, 'tall');
    orient(gcf, 'landscape');
end
end

%%
function [popCorrs, noiseCorrs] = ExaminePopAndNoiseCorrs(preTrial, rlyTrial, latTrial, pstTrial, sortVect, sortVectID)
preSorted = preTrial(sortVect,:);
rlySorted = rlyTrial(sortVect,:);
latSorted = latTrial(sortVect,:);
pstSorted = pstTrial(sortVect,:);

figure;
subplot(4,4,1)
noiseCorrs{1} = corr(preSorted);
imagesc(noiseCorrs{1}, [-1 1]);
title('Pre-Trial Period: Unit Noise Correlations');
subplot(4,4,2)
histogram(noiseCorrs{1}(logical(triu(ones(size(noiseCorrs{1})),1))), -1:0.01:1, 'edgecolor', 'none', 'facecolor', 'k');
set(gca, 'color', 'none', 'box', 'off');
grid on
hold on
% line([0 0], get(gca,'ylim'), 'linestyle', ':', 'color', 'r', 'linewidth', 1);
subplot(4,4,3)
popCorrs{1} = corr(preSorted');
imagesc(popCorrs{1}, [-1 1]);
title('Trial Population Correlations');
subplot(4,4,4)
histogram(popCorrs{1}(logical(triu(ones(size(popCorrs{1})),1))), -1:0.01:1, 'edgecolor', 'none', 'facecolor', 'k');
set(gca, 'color', 'none', 'box', 'off');
grid on
hold on
% line([0 0], get(gca,'ylim'), 'linestyle', ':', 'color', 'r', 'linewidth', 1);

subplot(4,4,5)
noiseCorrs{2} = corr(rlySorted);
imagesc(noiseCorrs{2}, [-1 1]);
title('Early-Trial Period: Unit Noise Correlations');
subplot(4,4,6)
histogram(noiseCorrs{2}(logical(triu(ones(size(noiseCorrs{2})),1))), -1:0.01:1, 'edgecolor', 'none', 'facecolor', 'k');
set(gca, 'color', 'none', 'box', 'off');
grid on
hold on
% line([0 0], get(gca,'ylim'), 'linestyle', ':', 'color', 'r', 'linewidth', 1);
subplot(4,4,7)
popCorrs{2} = corr(rlySorted');
imagesc(popCorrs{2}, [-1 1]);
title('Trial Population Correlations');
subplot(4,4,8)
histogram(popCorrs{2}(logical(triu(ones(size(popCorrs{2})),1))), -1:0.01:1, 'edgecolor', 'none', 'facecolor', 'k');
set(gca, 'color', 'none', 'box', 'off');
grid on
hold on
% line([0 0], get(gca,'ylim'), 'linestyle', ':', 'color', 'r', 'linewidth', 1);

subplot(4,4,9)
noiseCorrs{3} = corr(latSorted);
imagesc(noiseCorrs{3}, [-1 1]);
title('Late-Trial Period: Unit Noise Correlations');
subplot(4,4,10)
histogram(noiseCorrs{3}(logical(triu(ones(size(noiseCorrs{3})),1))), -1:0.01:1, 'edgecolor', 'none', 'facecolor', 'k');
set(gca, 'color', 'none', 'box', 'off');
grid on
hold on
% line([0 0], get(gca,'ylim'), 'linestyle', ':', 'color', 'r', 'linewidth', 1);
subplot(4,4,11)
popCorrs{3} = corr(latSorted');
imagesc(popCorrs{3}, [-1 1]);
title('Trial Population Correlations');
subplot(4,4,12)
histogram(popCorrs{3}(logical(triu(ones(size(popCorrs{3})),1))), -1:0.01:1, 'edgecolor', 'none', 'facecolor', 'k');
set(gca, 'color', 'none', 'box', 'off');
grid on
hold on
% line([0 0], get(gca,'ylim'), 'linestyle', ':', 'color', 'r', 'linewidth', 1);

subplot(4,4,13)
noiseCorrs{4} = corr(pstSorted);
imagesc(noiseCorrs{4}, [-1 1]);
title('Post-Trial Period: Unit Noise Correlations');
subplot(4,4,14)
histogram(noiseCorrs{4}(logical(triu(ones(size(noiseCorrs{4})),1))), -1:0.01:1, 'edgecolor', 'none', 'facecolor', 'k');
set(gca, 'color', 'none', 'box', 'off');
grid on
hold on
% line([0 0], get(gca,'ylim'), 'linestyle', ':', 'color', 'r', 'linewidth', 1);
subplot(4,4,15)
popCorrs{4} = corr(pstSorted');
imagesc(popCorrs{4}, [-1 1]);
title('Trial Population Correlations');
subplot(4,4,16)
histogram(popCorrs{4}(logical(triu(ones(size(popCorrs{4})),1))), -1:0.01:1, 'edgecolor', 'none', 'facecolor', 'k');
set(gca, 'color', 'none', 'box', 'off');
grid on
hold on
% line([0 0], get(gca,'ylim'), 'linestyle', ':', 'color', 'r', 'linewidth', 1);

annotation('textbox', 'position', [0.025 0.935 0.7 0.05], 'String', ['\bf\fontsize{14}' sortVectID],...
    'linestyle', 'none', 'horizontalalignment', 'left');
    
curDir = cd;
annotation('textbox', 'position', [0.025 0.025 0.7 0.05], 'String', curDir,...
    'linestyle', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');

colormap jet

axesHandles = findobj(get(gcf,'Children'), 'flat','Type','axes');
axis(axesHandles,'square')
orient(gcf, 'tall');
orient(gcf, 'landscape');
end


