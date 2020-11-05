%% PFCdualList_MLB
clc
clear all

%% Runtime variables
piWindow = [-0.5 0.5];
poWindow = [-0.5 0.5];
binSize = 200;
dsRate = 5;
cLim = [0 0.01];

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
% trialInfo = [1:Trial#, 2:Sequence#, 3:Position, 4:Odor, 5:Performance, 6:tDist, 7:PokeDur, 8:WithdrawLat]
[piSpkMtx, piTrialTime, trialInfo] = OrganizeAndBinSpikes(ensembleMatrix, behavMatrix, behavMatrixColIDs, 'PokeIn', piWindow, binSize, dsRate);
[poSpkMtx, poTrialTime, ~] = OrganizeAndBinSpikes(ensembleMatrix, behavMatrix, behavMatrixColIDs, 'PokeOut', poWindow, binSize, dsRate);

decodings.PItime = piTrialTime;
decodings.POtime = poTrialTime;

clear ensembleMatrix ensembleMatrixColIDs behavMatrix behavMatrixColIDs

%% Organize Data into Sequences & Identify Novel and Familiar Sequences
seqNums = unique(trialInfo(:,2));
seqCells = cell(size(seqNums));
seqTdist = zeros(size(seqNums));
for obsSeq = 1:length(seqNums)
    curSeq = trialInfo(trialInfo(:,2)==seqNums(obsSeq),:);
    seqCells{obsSeq} = curSeq;
    seqTdist(obsSeq) = sum(curSeq(:,6));
end
seqLength = cellfun(@(a)size(a,1),seqCells);
listID = cellfun(@(a)a(1,4),seqCells);

famFISs = seqCells(listID==1 & seqLength==4 & seqTdist==0)';
novFISs = seqCells(listID==11 & seqLength==4 & seqTdist==0)';

famFIStrls = cell2mat(cellfun(@(a)a(:,1), famFISs, 'uniformoutput', 0));
novFIStrls = cell2mat(cellfun(@(a)a(:,1), novFISs, 'uniformoutput', 0));

if ~isempty(intersect(famFIStrls, novFIStrls))
    error('Non-unique Familiar and Novel sequences, check code and files!');
end

isSeqs = [famFIStrls, novFIStrls];
issIDlog = [false(1,size(famFIStrls,2)), true(1,size(novFIStrls,2))];
%% Run decoding using leave-one-out
issPosts = nan((size(piSpkMtx,1) + size(poSpkMtx,1))*4, (size(piSpkMtx,1) + size(poSpkMtx,1))*8,size([famFIStrls, novFIStrls],2));
for s = 1:size(isSeqs,2)
    tempISSs = isSeqs;
    tempISSs(:,s) = [];
    tempISSidLog = issIDlog;
    tempISSidLog(s) = [];
    
    tempFAMfisLikes = cell(4,1);
    tempNOVfisLikes = cell(4,1);
    tempObsv = cell(4,1);
    for p = 1:4 
        tempFAMfisLikes{p} = [mean(piSpkMtx(:,:,tempISSs(p,~tempISSidLog)),3); mean(poSpkMtx(:,:,tempISSs(p,~tempISSidLog)),3)];
        tempNOVfisLikes{p} = [mean(piSpkMtx(:,:,tempISSs(p,tempISSidLog)),3); mean(poSpkMtx(:,:,tempISSs(p,tempISSidLog)),3)];
        tempObsv{p} = [piSpkMtx(:,:,isSeqs(p,s)); poSpkMtx(:,:,isSeqs(p,s))];
    end
    
    issPosts(:,:,s) = CalcStaticBayesPost(cell2mat([tempFAMfisLikes;tempNOVfisLikes]), cell2mat(tempObsv), binSize);
end

%% Plot LOO Data
famISSposts = issPosts(:,:,~issIDlog);
novISSposts = issPosts(:,:,issIDlog);
figure;
famSP = subplot(6,8,[1:4,9:12,17:20,25:28]);
imagesc(nanmean(famISSposts,3)', cLim); set(gca, 'ydir', 'normal'); colormap jet
set(gca, 'xtick', [], 'ytick', []);
for op = 1:8
    line([size(piSpkMtx,1) + size(poSpkMtx,1), size(piSpkMtx,1) + size(poSpkMtx,1)]*op, get(gca, 'ylim'), 'linestyle', '-', 'color', 'white', 'linewidth', 2);
    line([size(piSpkMtx,1)/2 size(piSpkMtx,1)/2] + repmat(size(piSpkMtx,1) + size(poSpkMtx,1),1,2)*(op-1), get(gca, 'ylim'), 'linestyle', ':', 'color', 'white', 'linewidth', 1);
    line([size(poSpkMtx,1)/2 + size(piSpkMtx,1)*op + size(poSpkMtx,1)*(op-1) size(poSpkMtx,1)/2 + size(piSpkMtx,1)*op + size(poSpkMtx,1)*(op-1)], get(gca, 'ylim'), 'linestyle', ':', 'color', 'white', 'linewidth', 1);
    line(get(gca, 'xlim'), [size(piSpkMtx,1) + size(poSpkMtx,1), size(piSpkMtx,1) + size(poSpkMtx,1)]*op, 'linestyle', '-', 'color', 'white', 'linewidth', 2);
    line(get(gca, 'xlim'), [size(piSpkMtx,1)/2 size(piSpkMtx,1)/2] + repmat(size(piSpkMtx,1) + size(poSpkMtx,1),1,2)*(op-1), 'linestyle', ':', 'color', 'white', 'linewidth', 1);
    line(get(gca, 'xlim'), [size(poSpkMtx,1)/2 + size(piSpkMtx,1)*op + size(poSpkMtx,1)*(op-1) size(poSpkMtx,1)/2 + size(piSpkMtx,1)*op + size(poSpkMtx,1)*(op-1)], 'linestyle', ':', 'color', 'white', 'linewidth', 1);
end
line(get(gca,'xlim'), [size(piSpkMtx,1)+size(poSpkMtx,1), size(piSpkMtx,1) + size(poSpkMtx,1)]*4, 'linestyle', '--', 'color', 'k', 'linewidth', 2);

novSP = subplot(6,8,[5:8,13:16,21:24,29:32]);
imagesc(nanmean(novISSposts,3)', cLim); set(gca, 'ydir', 'normal'); colormap jet
set(gca, 'xtick', [], 'ytick', []);
for op = 1:8
    line([size(piSpkMtx,1) + size(poSpkMtx,1), size(piSpkMtx,1) + size(poSpkMtx,1)]*op, get(gca, 'ylim'), 'linestyle', '-', 'color', 'white', 'linewidth', 2);
    line([size(piSpkMtx,1)/2 size(piSpkMtx,1)/2] + repmat(size(piSpkMtx,1) + size(poSpkMtx,1),1,2)*(op-1), get(gca, 'ylim'), 'linestyle', ':', 'color', 'white', 'linewidth', 1);
    line([size(poSpkMtx,1)/2 + size(piSpkMtx,1)*op + size(poSpkMtx,1)*(op-1) size(poSpkMtx,1)/2 + size(piSpkMtx,1)*op + size(poSpkMtx,1)*(op-1)], get(gca, 'ylim'), 'linestyle', ':', 'color', 'white', 'linewidth', 1);
    line(get(gca, 'xlim'), [size(piSpkMtx,1) + size(poSpkMtx,1), size(piSpkMtx,1) + size(poSpkMtx,1)]*op, 'linestyle', '-', 'color', 'white', 'linewidth', 2);
    line(get(gca, 'xlim'), [size(piSpkMtx,1)/2 size(piSpkMtx,1)/2] + repmat(size(piSpkMtx,1) + size(poSpkMtx,1),1,2)*(op-1), 'linestyle', ':', 'color', 'white', 'linewidth', 1);
    line(get(gca, 'xlim'), [size(poSpkMtx,1)/2 + size(piSpkMtx,1)*op + size(poSpkMtx,1)*(op-1) size(poSpkMtx,1)/2 + size(piSpkMtx,1)*op + size(poSpkMtx,1)*(op-1)], 'linestyle', ':', 'color', 'white', 'linewidth', 1);
end
line(get(gca,'xlim'), [size(piSpkMtx,1)+size(poSpkMtx,1), size(piSpkMtx,1) + size(poSpkMtx,1)]*4, 'linestyle', '--', 'color', 'k', 'linewidth', 2);

%% Decode LOO data
trialTemplate = ones(1,size(piSpkMtx,1)+size(poSpkMtx,1));
odorLog = [trialTemplate*1,...
    trialTemplate*2,...
    trialTemplate*3,...
    trialTemplate*4,...
    trialTemplate*11,...
    trialTemplate*12,...
    trialTemplate*13,...
    trialTemplate*14];
posLog = [trialTemplate*1,...
    trialTemplate*2,...
    trialTemplate*3,...
    trialTemplate*4,...
    trialTemplate*1,...
    trialTemplate*2,...
    trialTemplate*3,...
    trialTemplate*4];
seqLog = [repmat(trialTemplate*1,1,4),repmat(trialTemplate*2,1,4)];
timeLog = repmat([piTrialTime', poTrialTime'+1+dsRate/1000], 1,8);

% Decode Odor *************************************************************
odorDecode = DecodeBayesPost(issPosts, odorLog);

famOdorDecode = odorDecode(:,~issIDlog);
decodings.FamOdor = [mean(famOdorDecode==1,2), mean(famOdorDecode==2,2), mean(famOdorDecode==3,2), mean(famOdorDecode==4,2),...
    mean(famOdorDecode==11,2), mean(famOdorDecode==12,2), mean(famOdorDecode==13,2), mean(famOdorDecode==14,2)];
famOdrSP = subplot(6,8,33:36);
plot(1:size(odorDecode,1), smooth(decodings.FamOdor(:,1)-decodings.FamOdor(:,5)), 'color', [44/255 168/255 224/255], 'linewidth', 1);
hold on;
plot(1:size(odorDecode,1), smooth(decodings.FamOdor(:,2)-decodings.FamOdor(:,6)), 'color', [154/255 133/255 122/255], 'linewidth', 1);
plot(1:size(odorDecode,1), smooth(decodings.FamOdor(:,3)-decodings.FamOdor(:,7)), 'color', [9/255 161/255 74/255], 'linewidth', 1);
plot(1:size(odorDecode,1), smooth(decodings.FamOdor(:,4)-decodings.FamOdor(:,8)), 'color', [128/255 66/255 151/255], 'linewidth', 1);
legend('A-W', 'B-X', 'C-Y','D-Z', 'location', 'southoutside', 'orientation', 'horizontal');
ylabel([{'List Difference'};{'(Fam-Nov)'}]);
axis tight
set(gca, 'ylim', [-1 1], 'xtick', []);
for op = 1:4
    line([size(piSpkMtx,1) + size(poSpkMtx,1), size(piSpkMtx,1) + size(poSpkMtx,1)]*op, get(gca, 'ylim'), 'linestyle', '-', 'color', 'k', 'linewidth', 2);  
    line([size(piSpkMtx,1)/2 size(piSpkMtx,1)/2]+(size(piSpkMtx,1)+size(poSpkMtx,1))*(op-1), get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
    line([size(poSpkMtx,1)/2+size(piSpkMtx,1) size(poSpkMtx,1)/2+size(piSpkMtx,1)]+(size(piSpkMtx,1)+size(poSpkMtx,1))*(op-1), get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
end

novOdorDecode = odorDecode(:,issIDlog);
decodings.NovOdor = [mean(novOdorDecode==1,2), mean(novOdorDecode==2,2), mean(novOdorDecode==3,2), mean(novOdorDecode==4,2),...
    mean(novOdorDecode==11,2), mean(novOdorDecode==12,2), mean(novOdorDecode==13,2), mean(novOdorDecode==14,2)];
novOdrSP = subplot(6,8,37:40);
plot(1:size(odorDecode,1), smooth(decodings.NovOdor(:,1)-decodings.NovOdor(:,5)), 'color', [44/255 168/255 224/255], 'linewidth', 1);
hold on;
plot(1:size(odorDecode,1), smooth(decodings.NovOdor(:,2)-decodings.NovOdor(:,6)), 'color', [154/255 133/255 122/255], 'linewidth', 1);
plot(1:size(odorDecode,1), smooth(decodings.NovOdor(:,3)-decodings.NovOdor(:,7)), 'color', [9/255 161/255 74/255], 'linewidth', 1);
plot(1:size(odorDecode,1), smooth(decodings.NovOdor(:,4)-decodings.NovOdor(:,8)), 'color', [128/255 66/255 151/255], 'linewidth', 1);
legend('A-W', 'B-X', 'C-Y','D-Z', 'location', 'southoutside', 'orientation', 'horizontal');
axis tight
set(gca, 'ylim', [-1 1], 'xtick', []);
for op = 1:4
    line([size(piSpkMtx,1) + size(poSpkMtx,1), size(piSpkMtx,1) + size(poSpkMtx,1)]*op, get(gca, 'ylim'), 'linestyle', '-', 'color', 'k', 'linewidth', 2);  
    line([size(piSpkMtx,1)/2 size(piSpkMtx,1)/2]+(size(piSpkMtx,1)+size(poSpkMtx,1))*(op-1), get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
    line([size(poSpkMtx,1)/2+size(piSpkMtx,1) size(poSpkMtx,1)/2+size(piSpkMtx,1)]+(size(piSpkMtx,1)+size(poSpkMtx,1))*(op-1), get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
end

% Decode Position *********************************************************
posDecode = DecodeBayesPost(issPosts, posLog);

famPosDecode = posDecode(:,~issIDlog);
decodings.FamPos = [mean(famPosDecode==1,2), mean(famPosDecode==2,2), mean(famPosDecode==3,2), mean(famPosDecode==4,2)];
famPosSP = subplot(6,8,41:44);
plot(1:size(odorDecode,1), smooth(decodings.FamPos(:,1)), 'color', [44/255 168/255 224/255], 'linewidth', 1);
hold on;
plot(1:size(odorDecode,1), smooth(decodings.FamPos(:,2)), 'color', [154/255 133/255 122/255], 'linewidth', 1);
plot(1:size(odorDecode,1), smooth(decodings.FamPos(:,3)), 'color', [9/255 161/255 74/255], 'linewidth', 1);
plot(1:size(odorDecode,1), smooth(decodings.FamPos(:,4)), 'color', [128/255 66/255 151/255], 'linewidth', 1);
legend('1', '2', '3','4', 'location', 'southoutside', 'orientation', 'horizontal');
ylabel([{'Decoding'};{'(% Trials)'}]);
axis tight
set(gca, 'ylim', [0 1], 'xtick', []);
for op = 1:4
    line([size(piSpkMtx,1) + size(poSpkMtx,1), size(piSpkMtx,1) + size(poSpkMtx,1)]*op, get(gca, 'ylim'), 'linestyle', '-', 'color', 'k', 'linewidth', 2);  
    line([size(piSpkMtx,1)/2 size(piSpkMtx,1)/2]+(size(piSpkMtx,1)+size(poSpkMtx,1))*(op-1), get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
    line([size(poSpkMtx,1)/2+size(piSpkMtx,1) size(poSpkMtx,1)/2+size(piSpkMtx,1)]+(size(piSpkMtx,1)+size(poSpkMtx,1))*(op-1), get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
end

novPosDecode = posDecode(:,issIDlog);
decodings.NovPos = [mean(novPosDecode==1,2), mean(novPosDecode==2,2), mean(novPosDecode==3,2), mean(novPosDecode==4,2)];
novPosSP = subplot(6,8,45:48);
plot(1:size(odorDecode,1), smooth(decodings.NovPos(:,1)), 'color', [44/255 168/255 224/255], 'linewidth', 1);
hold on;
plot(1:size(odorDecode,1), smooth(decodings.NovPos(:,2)), 'color', [154/255 133/255 122/255], 'linewidth', 1);
plot(1:size(odorDecode,1), smooth(decodings.NovPos(:,3)), 'color', [9/255 161/255 74/255], 'linewidth', 1);
plot(1:size(odorDecode,1), smooth(decodings.NovPos(:,4)), 'color', [128/255 66/255 151/255], 'linewidth', 1);
legend('1', '2', '3','4', 'location', 'southoutside', 'orientation', 'horizontal');
axis tight
set(gca, 'ylim', [0 1], 'xtick', []);
for op = 1:4
    line([size(piSpkMtx,1) + size(poSpkMtx,1), size(piSpkMtx,1) + size(poSpkMtx,1)]*op, get(gca, 'ylim'), 'linestyle', '-', 'color', 'k', 'linewidth', 2);  
    line([size(piSpkMtx,1)/2 size(piSpkMtx,1)/2]+(size(piSpkMtx,1)+size(poSpkMtx,1))*(op-1), get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
    line([size(poSpkMtx,1)/2+size(piSpkMtx,1) size(poSpkMtx,1)/2+size(piSpkMtx,1)]+(size(piSpkMtx,1)+size(poSpkMtx,1))*(op-1), get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
end

annotation('textbox', 'position', [0 0.935 1 0.05], 'String', ['\bf\fontsize{10}' sprintf('Bin = %i ms; Step = %i ms; PokeInWindow = [%ims +%ims], PokeOutWindow = [%ims +%ims]', binSize, dsRate, piWindow(1)*1000,piWindow(2)*1000,poWindow(1)*1000, poWindow(2)*1000)],...
    'linestyle', 'none', 'horizontalalignment', 'right');
curDir = cd;
annotation('textbox', 'position', [0.025 0.025 0.7 0.05], 'String', curDir,...
    'linestyle', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');

figure;
subplot(1,2,1)
novDiff = decodings.NovOdor(:,1:4)-decodings.NovOdor(:,5:8);
ndCounts = histcounts(novDiff(:), -1:0.1:1);
famDiff = decodings.FamOdor(:,1:4)-decodings.FamOdor(:,5:8);
fdCounts = histcounts(famDiff(:), -1:0.1:1);
fB = bar(-0.95:0.1:0.95, fdCounts,1);
fB.FaceAlpha = 0.5;
hold on;
nB = bar(-0.95:0.1:0.95, ndCounts, 1, 'r');
nB.FaceAlpha = 0.5;
legend('Seq1', 'Seq2');
xlabel([{'Decoding Difference'}; {'Seq1-Seq2'}]);
subplot(1,2,2)
fB = plot(-0.95:0.1:0.95, cumsum(fdCounts)./sum(fdCounts));
hold on;
nB = plot(-0.95:0.1:0.95, cumsum(ndCounts)./sum(ndCounts), 'color', 'r');
xlabel([{'Decoding Difference'}; {'Seq1-Seq2'}]);
ylabel('Cumulative Proportion');

annotation('textbox', 'position', [0 0.935 1 0.05], 'String', ['\bf\fontsize{10}' sprintf('Bin = %i ms; Step = %i ms; PokeInWindow = [%ims +%ims], PokeOutWindow = [%ims +%ims]', binSize, dsRate, piWindow(1)*1000,piWindow(2)*1000,poWindow(1)*1000, poWindow(2)*1000)],...
    'linestyle', 'none', 'horizontalalignment', 'right');
curDir = cd;
annotation('textbox', 'position', [0.025 0.025 0.7 0.05], 'String', curDir,...
    'linestyle', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');

% Decode Time *************************************************************
timeDecode = DecodeBayesPost(issPosts, timeLog);

% Decode Sequence *********************************************************
seqDecode = DecodeBayesPost(issPosts, seqLog);

%% Examine Difference in Decoding Using pdfs: Nov X Fam
famFamPost = nanmean(issPosts(:,seqLog==1,~issIDlog),3);
famNovPost = nanmean(issPosts(:,seqLog==2,~issIDlog),3);
novFamPost = nanmean(issPosts(:,seqLog==1,issIDlog),3);
novNovPost = nanmean(issPosts(:,seqLog==2,issIDlog),3);

figure; 
subplot(3,4,2)
corrScatPlot(famFamPost(:), famNovPost(:), 'Fam-Fam', 'Fam-Nov', 'rmvZro');
subplot(3,4,3)
corrScatPlot(famFamPost(:), novFamPost(:), 'Fam-Fam', 'Nov-Fam', 'rmvZro');
subplot(3,4,4)
corrScatPlot(famFamPost(:), novNovPost(:), 'Fam-Fam', 'Nov-Nov', 'rmvZro');
subplot(3,4,7)
corrScatPlot(famNovPost(:), novFamPost(:), 'Fam-Nov', 'Nov-Fam', 'rmvZro');
subplot(3,4,8)
corrScatPlot(famNovPost(:), novNovPost(:), 'Fam-Nov', 'Nov-Nov', 'rmvZro');
subplot(3,4,12)
corrScatPlot(novFamPost(:), novNovPost(:), 'Nov-Fam', 'Nov-Nov', 'rmvZro');

annotation('textbox', 'position', [0 0.935 1 0.05], 'String', ['\bf\fontsize{10}' sprintf('Bin = %i ms; Step = %i ms; PokeInWindow = [%ims +%ims], PokeOutWindow = [%ims +%ims]', binSize, dsRate, piWindow(1)*1000,piWindow(2)*1000,poWindow(1)*1000, poWindow(2)*1000)],...
    'linestyle', 'none', 'horizontalalignment', 'right');
curDir = cd;
annotation('textbox', 'position', [0.025 0.025 0.7 0.05], 'String', curDir,...
    'linestyle', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');

figure;
subplot(2,2,1)
famDiff = famFamPost - famNovPost;
imagesc(famDiff, [-0.01, 0.01]); set(gca, 'ydir', 'normal'); colormap jet
title([{'Familiar Template'};{'Decode Fam - Decode Nov'}]);
subplot(2,2,2)
novDiff = novFamPost - novNovPost;
imagesc(novDiff, [-0.01, 0.01]); set(gca, 'ydir', 'normal'); colormap jet
title([{'Novel Template'};{'Decode Fam - Decode Nov'}]);
subplot(2,2,3)
corrScatPlot(famDiff(:), novDiff(:), 'FamDiff', 'NovDiff', 'rmvZro');
title('Difference Similarity');

annotation('textbox', 'position', [0 0.935 1 0.05], 'String', ['\bf\fontsize{10}' sprintf('Bin = %i ms; Step = %i ms; PokeInWindow = [%ims +%ims], PokeOutWindow = [%ims +%ims]', binSize, dsRate, piWindow(1)*1000,piWindow(2)*1000,poWindow(1)*1000, poWindow(2)*1000)],...
    'linestyle', 'none', 'horizontalalignment', 'right');
curDir = cd;
annotation('textbox', 'position', [0.025 0.025 0.7 0.05], 'String', curDir,...
    'linestyle', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');


%% Collapse Across Odor
odorCollapseDecode = zeros(4,8,2);
odors = unique(odorLog);
osbPosLog = [trialTemplate*1,...
    trialTemplate*2,...
    trialTemplate*3,...
    trialTemplate*4];
for obsSeq = 1:2
    curPosts = nanmean(issPosts(:,:,issIDlog==(obsSeq-1)),3);
    for obsPos = 1:4
        curObsPosLog = osbPosLog==obsPos;
        for decOdr = 1:length(odors)
            curDecOdrLog = odorLog==odors(decOdr);
            odorCollapseDecode(obsPos,decOdr,obsSeq) = mean(sum(curPosts(curObsPosLog,curDecOdrLog),2));
        end
    end
end

figure;
ffSP = subplot(3,3,1);
imagesc(odorCollapseDecode(:,1:4,1), [0 0.5]);
set(gca, 'xtick', 1:4, 'ytick', 1:4, 'xticklabel', [{'A'}, {'B'}, {'C'}, {'D'}], 'yticklabel', [{'A'}, {'B'}, {'C'}, {'D'}]);
colormap(ffSP, 'jet');
title('Fam Decode Fam');
fnSP = subplot(3,3,2);
imagesc(odorCollapseDecode(:,5:8,1), [0 0.5]);
set(gca, 'xtick', 1:4, 'ytick', 1:4, 'xticklabel', [{'A'}, {'B'}, {'C'}, {'D'}], 'yticklabel', [{'W'}, {'X'}, {'Y'}, {'X'}]);
colormap(fnSP, 'jet');
title('Fam Decode Nov');
fDiffSP = subplot(3,3,3);
famDecodeDiff = odorCollapseDecode(:,1:4,1) - odorCollapseDecode(:,5:8,1);
imagesc(famDecodeDiff, [-0.1 0.1]);
set(gca, 'xtick', 1:4, 'ytick', 1:4, 'xticklabel', [{'A-W'}, {'B-X'}, {'C-Y'}, {'D-Z'}], 'yticklabel', [{'W-A'}, {'X-B'}, {'Y-C'}, {'Z-D'}]);
colormap(fDiffSP, 'gray');
title('Fam-Nov Difference')
nfSP = subplot(3,3,4);
imagesc(odorCollapseDecode(:,1:4,2), [0 0.5]);
set(gca, 'xtick', 1:4, 'ytick', 1:4, 'xticklabel', [{'W'}, {'X'}, {'Y'}, {'X'}], 'yticklabel', [{'A'}, {'B'}, {'C'}, {'D'}]);
colormap(nfSP, 'jet');
title('Nov Decode Fam');
nnSP = subplot(3,3,5);
imagesc(odorCollapseDecode(:,5:8,2), [0 0.5]);
set(gca, 'xtick', 1:4, 'ytick', 1:4, 'xticklabel', [{'W'}, {'X'}, {'Y'}, {'X'}], 'yticklabel', [{'W'}, {'X'}, {'Y'}, {'X'}]);
colormap(nnSP, 'jet');
title('Nov Decode Nov');
nDiffSP = subplot(3,3,6);
novDecodeDiff = odorCollapseDecode(:,1:4,2)-odorCollapseDecode(:,5:8,2);
imagesc(novDecodeDiff, [-0.1 0.1]);
set(gca, 'xtick', 1:4, 'ytick', 1:4, 'xticklabel', [{'W-A'}, {'X-B'}, {'Y-C'}, {'Z-D'}], 'yticklabel', [{'A-W'}, {'B-X'}, {'C-Y'}, {'D-Z'}]);
colormap(nDiffSP, 'gray');

subplot(3,3,7)
tempFD = famDecodeDiff(2:4,2:4);
tempND = novDecodeDiff(2:4,2:4);
corrScatPlot(tempFD(:), tempND(:), 'Familiar Difference', 'Novel Difference');
title('Correlation Sans Pos1');
subplot(3,3,8)
tempFD = famDecodeDiff(:);
tempND = novDecodeDiff(:);
corrScatPlot(tempFD(2:end), tempND(2:end), 'Familiar Difference', 'Novel Difference');
title('Correlation Sans 1-1');
subplot(3,3,9)
corrScatPlot(tempFD(:), tempND(:), 'Familiar Difference', 'Novel Difference');
title('Correlation w/All');

            
%%
%%*******************************************************%%
%%********************** Functions **********************%%
%%*******************************************************%%
function [binnedMtx, trialTime, trialInfo] = OrganizeAndBinSpikes(ensembleMatrix, behavMatrix, behavMatrixColIDs, alignment, window, binSize, dsRate)
% First check to make sure reward signal is present within the
% behaviorMatrix
if sum(strcmp(behavMatrixColIDs, 'RewardSignal'))==0
    warning('Reward signal missing from behavior matrix, update the behavior matrix. OR ignore this if reward signal wasn''t recorded as plx flag');
end
% Organize behavior data and Extract Ensemble Spiking
% Taking 1/2 the binSize on either end to get rid of edge effects.
td = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [window(1)-(binSize/2/1000) window(2)+(binSize/2/1000)], alignment);
trialEnsemble = ExtractTrialData_SM(td, ensembleMatrix(:,2:end)); %#ok<*NODEF>
ensembleMtx = cell2mat(reshape(trialEnsemble, [1 1 length(trialEnsemble)]));
if strcmp(alignment, 'PokeIn')
    trialTimes = behavMatrix(td(1).TrialLogVect,1) - behavMatrix(td(1).PokeInIndex,1);
else
    trialTimes = behavMatrix(td(1).TrialLogVect,1) - behavMatrix(td(1).PokeOutIndex,1);
end

%% Bin the spiking data
% First convolve the entire trialEnsembleMtx with a square to bin the
% spikes
binnedEnsembleMtx = nan(size(ensembleMtx));
for t = 1:size(ensembleMtx,3)
    for u = 1:size(ensembleMtx,2)
        binnedEnsembleMtx(:,u,t) = conv(ensembleMtx(:,u,t), ones(1,binSize)./(binSize/1000), 'same');   % Spike Rate
%         binnedEnsembleMtx(:,u,t) = conv(trialEnsembleMtx(:,u,t), ones(1,binSize), 'same');    % Spike Count
    end
end
% Now remove the binSize/2 padding
unPaddedBinnedEnsembleMtx = binnedEnsembleMtx((binSize/2)+1:end-(binSize/2),:,:);
trialTimes = trialTimes((binSize/2)+1:end-(binSize/2));
% Now downsample the binned matrix
dsVect = downsample(1:size(unPaddedBinnedEnsembleMtx,1), dsRate);
binnedMtx = unPaddedBinnedEnsembleMtx(dsVect,:,:);
trialTime = trialTimes(dsVect);
   
trialInfo = [[td.TrialNum]', [td.SequenceNum]', [td.Position]', [td.Odor]', ...
    [td.Performance]',[td.TranspositionDistance]', [td.PokeDuration]', [td.WithdrawLatency]'];

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
function decode = TabulateBayesPost(post, id)
idS = unique(id);
decode = nan(size(post,1), size(post,3), length(idS));
for trl = 1:size(post,3)
    for t = 1:size(post,1)
        for iD = 1:length(idS)
            decode(t,trl,iD) = sum(post(t,id==idS(iD),trl));
        end
    end
end
end