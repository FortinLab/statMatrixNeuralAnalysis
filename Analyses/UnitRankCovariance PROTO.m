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

load(fileNames{cellfun(@(a)~isempty(a), strfind(fileNames, 'BehaviorMatrix'))});
load(fileNames{cellfun(@(a)~isempty(a), strfind(fileNames, 'EnsembleMatrix'))});

%% Extract Spike Data Based on Trial Period %%
pehBinSize = 0.01;

pokeInAlignedBehavMatrix = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [0 0], 'PokeIn');
pokeOutAlignedBehavMatrix = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [0 0], 'PokeOut');
trialBehavMatrix = CombineTrialLogVects_SM(pokeInAlignedBehavMatrix, pokeOutAlignedBehavMatrix);
trialEnsembles = ExtractTrialData_SM(trialBehavMatrix, ensembleMatrix(:,2:end));

binnedEnsembles = cell(size(trialEnsembles));
for trl = 1:length(trialEnsembles)
    [binnedEnsembles{trl}, pokeInBins] = RebinPEH2_SM(trialEnsembles{trl}, [], pehBinSize, 1000);
end

%% Calculate Unit Covariance Matrices & Assign Rank Values (Trial-Wise)
uniCovMat = cell(fliplr(size(binnedEnsembles)));
covMatUprLog = logical(triu(ones(size(ensembleMatrix,2)-1),1));
covMatNdx = reshape(1:(size(ensembleMatrix,2)-1)^2, [size(ensembleMatrix,2)-1 size(ensembleMatrix,2)-1]);
curMatUprNdx = covMatNdx(covMatUprLog);

for trl = 1:length(binnedEnsembles)
    tempCovMat = corr(binnedEnsembles{trl}');
    sortedTempCovMat = [tempCovMat(logical(triu(ones(size(tempCovMat)),1))), curMatUprNdx, nan(length(curMatUprNdx),1)];
    unqCovVals = unique(sortedTempCovMat(~isnan(sortedTempCovMat(:,1)),1));
    posCovVals = flipud(unqCovVals(unqCovVals>0));
    for pcv = 1:length(posCovVals)
        sortedTempCovMat(sortedTempCovMat(:,1)==posCovVals(pcv),3) = pcv;
    end
    negCovVals = unqCovVals(unqCovVals<0);
    for ncv = 1:length(negCovVals)
        sortedTempCovMat(sortedTempCovMat(:,1)==negCovVals(ncv),3) = ncv*-1;
    end
    uniCovMat{trl} = sortrows(sortedTempCovMat,-3);
end

% Evaluate Mean Rank Values
unpackedUCM = cell2mat(uniCovMat);
unpackedUCM(isnan(unpackedUCM(:,3)),:) = [];
reMtxdUCM = nan(size(covMatNdx));
reMtxdUCMvar = nan(size(covMatNdx));
reMtxdUCMpos = nan(size(covMatNdx));
reMtxdUCMposvar = nan(size(covMatNdx));
reMtxdUCMneg = nan(size(covMatNdx));
reMtxdUCMnegvar = nan(size(covMatNdx));
ndxs = unique(unpackedUCM(:,2));
for ndx = 1:length(unique(unpackedUCM(:,2)))
    curNdxDta = unpackedUCM(unpackedUCM(:,2)==ndxs(ndx),3);
    reMtxdUCM(ndxs(ndx)) = mean(curNdxDta);
    reMtxdUCMvar(ndxs(ndx)) = std(curNdxDta);
    reMtxdUCMpos(ndxs(ndx)) = mean(curNdxDta(curNdxDta>0));
    reMtxdUCMposvar(ndxs(ndx)) = std(curNdxDta(curNdxDta>0));
    reMtxdUCMneg(ndxs(ndx)) = mean(curNdxDta(curNdxDta<0))*-1;
    reMtxdUCMnegvar(ndxs(ndx)) = std(curNdxDta(curNdxDta<0));
end

% Plot Overall Vals
figure; 
subplot(4,4,1)
imagesc(reMtxdUCM);
hold on;
[r,c] = ind2sub(size(reMtxdUCM), find(isnan(reMtxdUCM)));
for p = 1:length(r)
    patch([c(p)-0.5 c(p)-0.5 c(p)+0.5 c(p)+0.5], [r(p)-0.5 r(p)+0.5 r(p)+0.5 r(p)-0.5], 'white', 'edgecolor', 'white');
end
set(gca, 'DataAspectRatio',[1 1 1],'Layer','top');
colorbar
title({'Mean Trial-Wise Covariance Rank Value' 'Closer to 0 = Higher Rank'});

subplot(4,4,2)
imagesc(reMtxdUCMpos);
hold on;
[r,c] = ind2sub(size(reMtxdUCMpos), find(isnan(reMtxdUCMpos)));
for p = 1:length(r)
    patch([c(p)-0.5 c(p)-0.5 c(p)+0.5 c(p)+0.5], [r(p)-0.5 r(p)+0.5 r(p)+0.5 r(p)-0.5], 'white', 'edgecolor', 'white');
end
set(gca, 'DataAspectRatio',[1 1 1],'Layer','top');
colorbar
title({'Mean Trial-Wise Covariance Rank Value' 'Positive Correlations Only'});

subplot(4,4,3)
imagesc(reMtxdUCMneg);
hold on;
[r,c] = ind2sub(size(reMtxdUCMneg), find(isnan(reMtxdUCMneg)));
for p = 1:length(r)
    patch([c(p)-0.5 c(p)-0.5 c(p)+0.5 c(p)+0.5], [r(p)-0.5 r(p)+0.5 r(p)+0.5 r(p)-0.5], 'white', 'edgecolor', 'white');
end
set(gca, 'DataAspectRatio',[1 1 1],'Layer','top');
colorbar
colormap jet
title({'Mean Trial-Wise Covariance Rank Value' 'Negative Correlations Only'});

% Plot Correlation between Positive and Negative Mean Vals
subplot(4,4,4)
corrScatPlot(reMtxdUCMpos((~isnan(reMtxdUCMneg) & ~isnan(reMtxdUCMpos))), reMtxdUCMneg((~isnan(reMtxdUCMneg) & ~isnan(reMtxdUCMpos))), 'Positive r-Val Trials Ranked', 'Negative r-Val Trials Ranked', 'Relationship of Pos and Neg r-Val Ranks');
set(gca, 'xlim', [min([get(gca, 'ylim') get(gca, 'xlim')]), max([get(gca, 'ylim') get(gca, 'xlim')])],...
    'ylim', [min([get(gca, 'ylim') get(gca, 'xlim')]), max([get(gca, 'ylim') get(gca, 'xlim')])],...
    'DataAspectRatio',[1 1 1],'Layer','top');

subplot(4,4,5)
imagesc(reMtxdUCMvar);
hold on;
[r,c] = ind2sub(size(reMtxdUCMvar), find(isnan(reMtxdUCMvar)));
for p = 1:length(r)
    patch([c(p)-0.5 c(p)-0.5 c(p)+0.5 c(p)+0.5], [r(p)-0.5 r(p)+0.5 r(p)+0.5 r(p)-0.5], 'white', 'edgecolor', 'white');
end
set(gca, 'DataAspectRatio',[1 1 1],'Layer','top');
colorbar
title({'Variance of Trial-Wise r-Val Rank Value' 'Closer to 0 = Higher Rank'});

subplot(4,4,6)
imagesc(reMtxdUCMposvar);
hold on;
[r,c] = ind2sub(size(reMtxdUCMposvar), find(isnan(reMtxdUCMposvar)));
for p = 1:length(r)
    patch([c(p)-0.5 c(p)-0.5 c(p)+0.5 c(p)+0.5], [r(p)-0.5 r(p)+0.5 r(p)+0.5 r(p)-0.5], 'white', 'edgecolor', 'white');
end
set(gca, 'DataAspectRatio',[1 1 1],'Layer','top');
colorbar
title({'Variance of Trial-Wise r-Val Rank Value' 'Positive Correlations Only'});

subplot(4,4,7)
imagesc(reMtxdUCMnegvar);
hold on;
[r,c] = ind2sub(size(reMtxdUCMnegvar), find(isnan(reMtxdUCMnegvar)));
for p = 1:length(r)
    patch([c(p)-0.5 c(p)-0.5 c(p)+0.5 c(p)+0.5], [r(p)-0.5 r(p)+0.5 r(p)+0.5 r(p)-0.5], 'white', 'edgecolor', 'white');
end
set(gca, 'DataAspectRatio',[1 1 1],'Layer','top');
colorbar
colormap jet
title({'Variance of Trial-Wise r-Val Rank Value' 'Negative Correlations Only'});

% Plot Correlation between Positive and Negative Vals Variance
subplot(4,4,8)
corrScatPlot(reMtxdUCMposvar((~isnan(reMtxdUCMnegvar) & ~isnan(reMtxdUCMposvar))), reMtxdUCMnegvar((~isnan(reMtxdUCMnegvar) & ~isnan(reMtxdUCMposvar))),...
    'Positive r-Val Rank Variance', 'Negative r-Val Rank Variance', 'Relationship of Variance of Pos and Neg r-Val Ranks');
set(gca, 'xlim', [min([get(gca, 'ylim') get(gca, 'xlim')]), max([get(gca, 'ylim') get(gca, 'xlim')])],...
    'ylim', [min([get(gca, 'ylim') get(gca, 'xlim')]), max([get(gca, 'ylim') get(gca, 'xlim')])],...
    'DataAspectRatio',[1 1 1],'Layer','top');

% Plot Correlation between Mean & Std Vals
subplot(4,4,9)
corrScatPlot(reMtxdUCM(~isnan(reMtxdUCM)), reMtxdUCMvar(~isnan(reMtxdUCM)), 'Rank Mean', 'Rank Variance', {'Mean vs Variance' 'Average Rank'});
set(gca, 'xlim', [min([get(gca, 'ylim') get(gca, 'xlim')]), max([get(gca, 'ylim') get(gca, 'xlim')])],...
    'ylim', [min([get(gca, 'ylim') get(gca, 'xlim')]), max([get(gca, 'ylim') get(gca, 'xlim')])],...
    'DataAspectRatio',[1 1 1],'Layer','top');

subplot(4,4,10)
corrScatPlot(reMtxdUCMpos(~isnan(reMtxdUCMpos)), reMtxdUCMposvar(~isnan(reMtxdUCMpos)), 'Rank Mean', 'Rank Variance', {'Mean vs Variance' 'Positive Values Rank'});
set(gca, 'xlim', [min([get(gca, 'ylim') get(gca, 'xlim')]), max([get(gca, 'ylim') get(gca, 'xlim')])],...
    'ylim', [min([get(gca, 'ylim') get(gca, 'xlim')]), max([get(gca, 'ylim') get(gca, 'xlim')])],...
    'DataAspectRatio',[1 1 1],'Layer','top');

subplot(4,4,11)
corrScatPlot(reMtxdUCMneg(~isnan(reMtxdUCMneg)), reMtxdUCMnegvar(~isnan(reMtxdUCMneg)), 'Rank Mean', 'Rank Variance', {'Mean vs Variance' 'Negative Values Rank'});
set(gca, 'xlim', [min([get(gca, 'ylim') get(gca, 'xlim')]), max([get(gca, 'ylim') get(gca, 'xlim')])],...
    'ylim', [min([get(gca, 'ylim') get(gca, 'xlim')]), max([get(gca, 'ylim') get(gca, 'xlim')])],...
    'DataAspectRatio',[1 1 1],'Layer','top');

subplot(4,4,13)
histogram(reMtxdUCM(~isnan(reMtxdUCM)));
title({'Distribution of Mean Overall Rank' '1 = Highest Rank'})
subplot(4,4,14)
histogram(reMtxdUCMpos(~isnan(reMtxdUCMpos)));
title('Distribution of Mean Positive r-Val Trial Rank')
subplot(4,4,15)
histogram(reMtxdUCMneg(~isnan(reMtxdUCMneg)));
title('Distribution of Mean Negativer-Val Trial Rank')

%% Calculate Unit Covariance Matrices & Assign Rank Values (Session-Wise)
uniCovMat = cell(fliplr(size(binnedEnsembles)));
covMatUprLog = logical(triu(ones(size(ensembleMatrix,2)-1),1));
covMatNdx = reshape(1:(size(ensembleMatrix,2)-1)^2, [size(ensembleMatrix,2)-1 size(ensembleMatrix,2)-1]);
curMatUprNdx = covMatNdx(covMatUprLog);

tempCovMat = corr(cell2mat(binnedEnsembles)');
sortedTempCovMat = [tempCovMat(logical(triu(ones(size(tempCovMat)),1))), curMatUprNdx, nan(length(curMatUprNdx),1)];
unqCovVals = unique(sortedTempCovMat(~isnan(sortedTempCovMat(:,1)),1));
posCovVals = flipud(unqCovVals(unqCovVals>0));
for pcv = 1:length(posCovVals)
    sortedTempCovMat(sortedTempCovMat(:,1)==posCovVals(pcv),3) = pcv;
end
negCovVals = unqCovVals(unqCovVals<0);
for ncv = 1:length(negCovVals)
    sortedTempCovMat(sortedTempCovMat(:,1)==negCovVals(ncv),3) = ncv*-1;
end
uniCovMat = sortrows(sortedTempCovMat,-3);

% Evaluate Mean Rank Values
unpackedUCM = uniCovMat;
unpackedUCM(isnan(unpackedUCM(:,3)),:) = [];
reMtxdUCM = nan(size(covMatNdx));
ndxs = unique(unpackedUCM(:,2));
for ndx = 1:length(unique(unpackedUCM(:,2)))
    curNdxDta = unpackedUCM(unpackedUCM(:,2)==ndxs(ndx),3);
    reMtxdUCM(ndxs(ndx)) = mean(curNdxDta);
end

% Plot CovVals
figure;
subplot(2,2,[1 3])
imagesc(reMtxdUCM, [max(abs(reMtxdUCM(:)))*-1 max(abs(reMtxdUCM(:)))]);
hold on;
[r,c] = ind2sub(size(reMtxdUCM), find(isnan(reMtxdUCM)));
for p = 1:length(r)
    patch([c(p)-0.5 c(p)-0.5 c(p)+0.5 c(p)+0.5], [r(p)-0.5 r(p)+0.5 r(p)+0.5 r(p)-0.5], 'white', 'edgecolor', 'white');
end
set(gca, 'DataAspectRatio',[1 1 1],'Layer','top');
colorbar
colormap jet
title('Session-Wide Unit Covariance Rank Values');

subplot(2,2,2)
histogram(reMtxdUCM(~isnan(reMtxdUCM)),[max(abs(reMtxdUCM(:)))*-1:5:max(abs(reMtxdUCM(:)))]);
set(gca, 'xlim', [max(abs(reMtxdUCM(:)))*-1 max(abs(reMtxdUCM(:)))]);
title('Rank Order Value Distribution')

subplot(2,2,4)
histogram(uniCovMat(:,1))
title('Overall Correlation Value Distribution')

