%%
clc
clear all

%% Runtime variables
piWindow = [-0.5 0.5];
poWindow = [-0.5 0.5];
binSize = 200;
dsRate = 5;
cLims = [0 0.005];

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
% trialInfo = [Trial#, Sequence#, Position, Odor, Performance, tDist, PokeDur, WithdrawLat]
[piSpkMtx, piTrialTime, trialInfo] = OrganizeAndBinSpikes(ensembleMatrix, behavMatrix, behavMatrixColIDs, 'PokeIn', piWindow, binSize, dsRate);
[poSpkMtx, poTrialTime, ~] = OrganizeAndBinSpikes(ensembleMatrix, behavMatrix, behavMatrixColIDs, 'PokeOut', poWindow, binSize, dsRate);

clear ensembleMatrix ensembleMatrixColIDs behavMatrix behavMatrixColIDs
%% Identify Fully InSeq Trials & Calculate FIS Liklihoods (PSTH)
fullInSeqSeqsStart = find(conv(trialInfo(:,4), 1:4, 'valid')==20 & conv(trialInfo(:,3), 1:4, 'valid')==20 & conv(trialInfo(:,5), ones(1,4), 'valid')==4);
inSeqSeqs = nan(3,length(fullInSeqSeqsStart));
for iS = 1:length(fullInSeqSeqsStart)
    inSeqSeqs(1,iS) = fullInSeqSeqsStart(iS);
    inSeqSeqs(2,iS) = fullInSeqSeqsStart(iS) + 1;
    inSeqSeqs(3,iS) = fullInSeqSeqsStart(iS) + 2;
    inSeqSeqs(4,iS) = fullInSeqSeqsStart(iS) + 3;
end

piFISlikes = repmat({nan(size(piSpkMtx,1), size(piSpkMtx,2))}, 4,1);
poFISlikes = repmat({nan(size(poSpkMtx,1), size(poSpkMtx,2))}, 4,1);
for p = 1:4
    piFISlikes{p} = mean(piSpkMtx(:,:,inSeqSeqs(p,:)),3);
    poFISlikes{p} = mean(poSpkMtx(:,:,inSeqSeqs(p,:)),3);
end

% Plot Likelihoods
figure; 
figLims = [0 max(max([cell2mat(piFISlikes); cell2mat(poFISlikes)]))];
subplot(4,2,1); imagesc(piTrialTime, 1:size(piFISlikes{1},2), piFISlikes{1}', figLims); title(gca, 'Poke In'); ylabel('Odor A'); hold on; line([0 0], get(gca, 'ylim'), 'color', 'w');
subplot(4,2,2); imagesc(poTrialTime, 1:size(poFISlikes{1},2), poFISlikes{1}', figLims); title(gca, 'Poke Out'); hold on; line([0 0], get(gca, 'ylim'), 'color', 'w');
subplot(4,2,3); imagesc(piTrialTime, 1:size(piFISlikes{1},2), piFISlikes{2}', figLims); ylabel('Odor B'); hold on; line([0 0], get(gca, 'ylim'), 'color', 'w');
subplot(4,2,4); imagesc(poTrialTime, 1:size(poFISlikes{1},2), poFISlikes{2}', figLims);  hold on; line([0 0], get(gca, 'ylim'), 'color', 'w');
subplot(4,2,5); imagesc(piTrialTime, 1:size(piFISlikes{1},2), piFISlikes{3}', figLims); ylabel('Odor C'); hold on; line([0 0], get(gca, 'ylim'), 'color', 'w');
subplot(4,2,6); imagesc(poTrialTime, 1:size(poFISlikes{1},2), poFISlikes{3}', figLims);  hold on; line([0 0], get(gca, 'ylim'), 'color', 'w');
subplot(4,2,7); imagesc(piTrialTime, 1:size(piFISlikes{1},2), piFISlikes{4}', figLims); ylabel('Odor D'); hold on; line([0 0], get(gca, 'ylim'), 'color', 'w');
subplot(4,2,8); imagesc(poTrialTime, 1:size(poFISlikes{1},2), poFISlikes{4}', figLims); hold on; line([0 0], get(gca, 'ylim'), 'color', 'w');

annotation('textbox', 'position', [0 0.935 1 0.05], 'String', ['\bf\fontsize{10}' sprintf('Bin = %i ms; Step = %i ms; PokeInWindow = [%ims +%ims], PokeOutWindow = [%ims +%ims]', binSize, dsRate, piWindow(1)*1000,piWindow(2)*1000,poWindow(1)*1000, poWindow(2)*1000)],...
    'linestyle', 'none', 'horizontalalignment', 'right');
curDir = cd;
annotation('textbox', 'position', [0.025 0.025 0.7 0.05], 'String', curDir,...
    'linestyle', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');

%% Decode via Leave 1 out
issPosts = nan((size(piSpkMtx,1) + size(poSpkMtx,1))*4, (size(piSpkMtx,1) + size(poSpkMtx,1))*4,size(inSeqSeqs,2));
for s = 1:size(inSeqSeqs,2)
    tempISSs = inSeqSeqs;
    tempISSs(:,s) = [];
    
    tempFISlikes = repmat({nan(size(piSpkMtx,1) + size(poSpkMtx,1), size(piSpkMtx,2) + size(poSpkMtx,2))}, 4,1);
    tempObsv = tempFISlikes;
    for p = 1:4 
        tempFISlikes{p} = [mean(piSpkMtx(:,:,tempISSs(p,:)),3); mean(poSpkMtx(:,:,tempISSs(p,:)),3)];
        tempObsv{p} = [piSpkMtx(:,:,inSeqSeqs(p,s)); poSpkMtx(:,:,inSeqSeqs(p,s))];
    end
    
    issPosts(:,:,s) = CalcStaticBayesPost(cell2mat(tempFISlikes), cell2mat(tempObsv), binSize);
end
odorLog = [ones(1,size(piSpkMtx,1)+size(poSpkMtx,1))*1,...
    ones(1,size(piSpkMtx,1)+size(poSpkMtx,1))*2,...
    ones(1,size(piSpkMtx,1)+size(poSpkMtx,1))*3,...
    ones(1,size(piSpkMtx,1)+size(poSpkMtx,1))*4];
timeLog = [piTrialTime', poTrialTime'+1+dsRate/1000,...
    piTrialTime', poTrialTime'+1+dsRate/1000,...
    piTrialTime', poTrialTime'+1+dsRate/1000,...
    piTrialTime', poTrialTime'+1+dsRate/1000];
    
figure; 
subplot(6,4,1:16);
% imagesc(nanmedian(issPosts,3)', [0 0.01]); set(gca, 'ydir', 'normal'); colormap jet
imagesc(nanmean(issPosts,3)', [0 0.01]); set(gca, 'ydir', 'normal'); colormap jet
set(gca, 'xtick', [], 'ytick', []);
line([size(piSpkMtx,1) + size(poSpkMtx,1), size(piSpkMtx,1) + size(poSpkMtx,1)], get(gca, 'ylim'), 'linestyle', '-', 'color', 'white', 'linewidth', 2);   
line([size(piSpkMtx,1)/2 size(piSpkMtx,1)/2], get(gca, 'ylim'), 'linestyle', ':', 'color', 'white', 'linewidth', 1);
line([size(poSpkMtx,1)/2+size(piSpkMtx,1) size(poSpkMtx,1)/2+size(piSpkMtx,1)], get(gca, 'ylim'), 'linestyle', ':', 'color', 'white', 'linewidth', 1);
line(get(gca, 'xlim'), [size(piSpkMtx,1) + size(poSpkMtx,1), size(piSpkMtx,1) + size(poSpkMtx,1)], 'linestyle', '-', 'color', 'white', 'linewidth', 2); 
line(get(gca, 'xlim'), [size(piSpkMtx,1)/2 size(piSpkMtx,1)/2], 'linestyle', ':', 'color', 'white', 'linewidth', 1);
line(get(gca, 'xlim'), [size(poSpkMtx,1)/2+size(piSpkMtx,1) size(poSpkMtx,1)/2+size(piSpkMtx,1)], 'linestyle', ':', 'color', 'white', 'linewidth', 1);

line([size(piSpkMtx,1) + size(poSpkMtx,1), size(piSpkMtx,1) + size(poSpkMtx,1)]*2, get(gca, 'ylim'), 'linestyle', '-', 'color', 'white', 'linewidth', 2);    
line([size(piSpkMtx,1)/2+size(piSpkMtx,1)+size(poSpkMtx,1) size(piSpkMtx,1)/2+size(piSpkMtx,1)+size(poSpkMtx,1)], get(gca, 'ylim'), 'linestyle', ':', 'color', 'white', 'linewidth', 1);
line([size(poSpkMtx,1)/2+size(piSpkMtx,1)*2+size(poSpkMtx,1) size(poSpkMtx,1)/2+size(piSpkMtx,1)*2+size(poSpkMtx,1)], get(gca, 'ylim'), 'linestyle', ':', 'color', 'white', 'linewidth', 1);     
line(get(gca, 'xlim'), [size(piSpkMtx,1) + size(poSpkMtx,1), size(piSpkMtx,1) + size(poSpkMtx,1)]*2, 'linestyle', '-', 'color', 'white', 'linewidth', 2);  
line(get(gca, 'xlim'), [size(piSpkMtx,1)/2+size(piSpkMtx,1)+size(poSpkMtx,1) size(piSpkMtx,1)/2+size(piSpkMtx,1)+size(poSpkMtx,1)], 'linestyle', ':', 'color', 'white', 'linewidth', 1);
line(get(gca, 'xlim'), [size(poSpkMtx,1)/2+size(piSpkMtx,1)*2+size(poSpkMtx,1) size(poSpkMtx,1)/2+size(piSpkMtx,1)*2+size(poSpkMtx,1)], 'linestyle', ':', 'color', 'white', 'linewidth', 1);

line([size(piSpkMtx,1) + size(poSpkMtx,1), size(piSpkMtx,1) + size(poSpkMtx,1)]*3, get(gca, 'ylim'), 'linestyle', '-', 'color', 'white', 'linewidth', 2); 
line([size(piSpkMtx,1)/2+size(piSpkMtx,1)*2+size(poSpkMtx,1)*2 size(piSpkMtx,1)/2+size(piSpkMtx,1)*2+size(poSpkMtx,1)*2], get(gca, 'ylim'), 'linestyle', ':', 'color', 'white', 'linewidth', 1);
line([size(poSpkMtx,1)/2+size(piSpkMtx,1)*3+size(poSpkMtx,1)*2 size(poSpkMtx,1)/2+size(piSpkMtx,1)*3+size(poSpkMtx,1)*2], get(gca, 'ylim'), 'linestyle', ':', 'color', 'white', 'linewidth', 1);           
line(get(gca, 'xlim'), [size(piSpkMtx,1) + size(poSpkMtx,1), size(piSpkMtx,1) + size(poSpkMtx,1)]*3, 'linestyle', '-', 'color', 'white', 'linewidth', 2);  
line(get(gca, 'xlim'), [size(piSpkMtx,1)/2+size(piSpkMtx,1)*2+size(poSpkMtx,1)*2 size(piSpkMtx,1)/2+size(piSpkMtx,1)*2+size(poSpkMtx,1)*2], 'linestyle', ':', 'color', 'white', 'linewidth', 1);
line(get(gca, 'xlim'), [size(poSpkMtx,1)/2+size(piSpkMtx,1)*3+size(poSpkMtx,1)*2 size(poSpkMtx,1)/2+size(piSpkMtx,1)*3+size(poSpkMtx,1)*2], 'linestyle', ':', 'color', 'white', 'linewidth', 1);

line([size(piSpkMtx,1) + size(poSpkMtx,1), size(piSpkMtx,1) + size(poSpkMtx,1)]*4, get(gca, 'ylim'), 'linestyle', '-', 'color', 'white', 'linewidth', 2);   
line([size(piSpkMtx,1)/2+size(piSpkMtx,1)*3+size(poSpkMtx,1)*3 size(piSpkMtx,1)/2+size(piSpkMtx,1)*3+size(poSpkMtx,1)*3], get(gca, 'ylim'), 'linestyle', ':', 'color', 'white', 'linewidth', 1);
line([size(poSpkMtx,1)/2+size(piSpkMtx,1)*4+size(poSpkMtx,1)*3 size(poSpkMtx,1)/2+size(piSpkMtx,1)*4+size(poSpkMtx,1)*3], get(gca, 'ylim'), 'linestyle', ':', 'color', 'white', 'linewidth', 1);     
line(get(gca, 'xlim'), [size(piSpkMtx,1) + size(poSpkMtx,1), size(piSpkMtx,1) + size(poSpkMtx,1)]*4, 'linestyle', '-', 'color', 'white', 'linewidth', 2);        
line(get(gca, 'xlim'), [size(piSpkMtx,1)/2+size(piSpkMtx,1)*3+size(poSpkMtx,1)*3 size(piSpkMtx,1)/2+size(piSpkMtx,1)*3+size(poSpkMtx,1)*3], 'linestyle', ':', 'color', 'white', 'linewidth', 1);
line(get(gca, 'xlim'), [size(poSpkMtx,1)/2+size(piSpkMtx,1)*4+size(poSpkMtx,1)*3 size(poSpkMtx,1)/2+size(piSpkMtx,1)*4+size(poSpkMtx,1)*3], 'linestyle', ':', 'color', 'white', 'linewidth', 1);

% % Odor Decoding (average posteriors)
% post = TabulateBayesPost(issPosts, odorLog);
% subplot(6,4,17:20)
% plot(1:size(post,1), nanmean(post(:,:,1),2), 'color', [44/255 168/255 224/255], 'linewidth', 1);
% hold on;
% plot(1:size(post,1), nanmean(post(:,:,2),2), 'color', [154/255 133/255 122/255], 'linewidth', 1);
% plot(1:size(post,1), nanmean(post(:,:,3),2), 'color', [9/255 161/255 74/255], 'linewidth', 1);
% plot(1:size(post,1), nanmean(post(:,:,4),2), 'color', [128/255 66/255 151/255], 'linewidth', 1);
% legend('A', 'B', 'C','D', 'location', 'southoutside', 'orientation', 'horizontal');
% ylabel([{'Odor Posterior'};{'(Mean +/- SEM)'}]);
% for o = 1:4
%     tempPostMean = nanmean(post(:,:,o),2);
% %     tempPostVar = nanstd(post(:,:,o),1,2);
%     tempPostVar = SEMcalc(post(:,:,o)')';
%     tp = patch('YData', [tempPostMean+tempPostVar; flipud(tempPostMean-tempPostVar)],...
%         'XData', [1:length(tempPostMean), length(tempPostMean):-1:1], 'FaceAlpha', .3);
%     if o==1
%         set(tp, 'FaceColor', [44/255 168/255 224/255], 'edgecolor', [44/255 168/255 224/255], 'EdgeAlpha', .5);
%     elseif o==2
%         set(tp, 'FaceColor', [154/255 133/255 122/255], 'edgecolor', [154/255 133/255 122/255], 'EdgeAlpha', .5);
%     elseif o==3
%         set(tp, 'FaceColor', [9/255 161/255 74/255], 'edgecolor', [9/255 161/255 74/255], 'EdgeAlpha', .5);
%     else
%         set(tp, 'FaceColor', [128/255 66/255 151/255], 'edgecolor', [128/255 66/255 151/255], 'EdgeAlpha', .5);
%     end 
% end

% Odor Decoding (true decoding)
decodeOdor = DecodeBayesPost(issPosts, odorLog);
subplot(6,4,17:20)
plot(1:size(decodeOdor,1), mean(decodeOdor==1,2), 'color', [44/255 168/255 224/255], 'linewidth', 1);
hold on;
plot(1:size(decodeOdor,1), mean(decodeOdor==2,2), 'color', [154/255 133/255 122/255], 'linewidth', 1);
plot(1:size(decodeOdor,1), mean(decodeOdor==3,2), 'color', [9/255 161/255 74/255], 'linewidth', 1);
plot(1:size(decodeOdor,1), mean(decodeOdor==4,2), 'color', [128/255 66/255 151/255], 'linewidth', 1);
legend('A', 'B', 'C','D', 'location', 'southoutside', 'orientation', 'horizontal');
ylabel([{'Decoding'};{'(% Trials)'}]);

axis tight
set(gca, 'ylim', [0 1], 'xtick', []);
line([size(piSpkMtx,1) + size(poSpkMtx,1), size(piSpkMtx,1) + size(poSpkMtx,1)], get(gca, 'ylim'), 'linestyle', '-', 'color', 'k', 'linewidth', 2);   
line([size(piSpkMtx,1) + size(poSpkMtx,1), size(piSpkMtx,1) + size(poSpkMtx,1)]*2, get(gca, 'ylim'), 'linestyle', '-', 'color', 'k', 'linewidth', 2);       
line([size(piSpkMtx,1) + size(poSpkMtx,1), size(piSpkMtx,1) + size(poSpkMtx,1)]*3, get(gca, 'ylim'), 'linestyle', '-', 'color', 'k', 'linewidth', 2);       
line([size(piSpkMtx,1) + size(poSpkMtx,1), size(piSpkMtx,1) + size(poSpkMtx,1)]*4, get(gca, 'ylim'), 'linestyle', '-', 'color', 'k', 'linewidth', 2);       
line([size(piSpkMtx,1)/2 size(piSpkMtx,1)/2], get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
line([size(poSpkMtx,1)/2+size(piSpkMtx,1) size(poSpkMtx,1)/2+size(piSpkMtx,1)], get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
line([size(piSpkMtx,1)/2+size(piSpkMtx,1)+size(poSpkMtx,1) size(piSpkMtx,1)/2+size(piSpkMtx,1)+size(poSpkMtx,1)], get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
line([size(poSpkMtx,1)/2+size(piSpkMtx,1)*2+size(poSpkMtx,1) size(poSpkMtx,1)/2+size(piSpkMtx,1)*2+size(poSpkMtx,1)], get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);  
line([size(piSpkMtx,1)/2+size(piSpkMtx,1)*2+size(poSpkMtx,1)*2 size(piSpkMtx,1)/2+size(piSpkMtx,1)*2+size(poSpkMtx,1)*2], get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
line([size(poSpkMtx,1)/2+size(piSpkMtx,1)*3+size(poSpkMtx,1)*2 size(poSpkMtx,1)/2+size(piSpkMtx,1)*3+size(poSpkMtx,1)*2], get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);       
line([size(piSpkMtx,1)/2+size(piSpkMtx,1)*3+size(poSpkMtx,1)*3 size(piSpkMtx,1)/2+size(piSpkMtx,1)*3+size(poSpkMtx,1)*3], get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
line([size(poSpkMtx,1)/2+size(piSpkMtx,1)*4+size(poSpkMtx,1)*3 size(poSpkMtx,1)/2+size(piSpkMtx,1)*4+size(poSpkMtx,1)*3], get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);     

% Temporal Decoding
decodeTime = DecodeBayesPost(issPosts, timeLog);
decodeLag = nan(size(decodeTime));
for c = 1:size(decodeTime,2)
    decodeLag(:,c) = decodeTime(:,c)-timeLog';
end
lagMean = nanmean(decodeLag,2);
lagVar = nanstd(decodeLag,1,2);
subplot(6,4,21:24)
plot(lagMean, '-k', 'linewidth', 1);
patch('YData', [lagMean+lagVar; flipud(lagMean-lagVar)],...
    'XData', [1:length(lagMean), length(lagMean):-1:1], 'FaceColor', 'k', 'FaceAlpha', .3, 'edgecolor', 'none');
axis tight
line([0 length(lagMean)], [0 0], 'linestyle', '--', 'color', 'k', 'linewidth', 1);
line([size(piSpkMtx,1) + size(poSpkMtx,1), size(piSpkMtx,1) + size(poSpkMtx,1)], get(gca, 'ylim'), 'linestyle', '-', 'color', 'k', 'linewidth', 2);   
line([size(piSpkMtx,1) + size(poSpkMtx,1), size(piSpkMtx,1) + size(poSpkMtx,1)]*2, get(gca, 'ylim'), 'linestyle', '-', 'color', 'k', 'linewidth', 2);       
line([size(piSpkMtx,1) + size(poSpkMtx,1), size(piSpkMtx,1) + size(poSpkMtx,1)]*3, get(gca, 'ylim'), 'linestyle', '-', 'color', 'k', 'linewidth', 2);       
line([size(piSpkMtx,1) + size(poSpkMtx,1), size(piSpkMtx,1) + size(poSpkMtx,1)]*4, get(gca, 'ylim'), 'linestyle', '-', 'color', 'k', 'linewidth', 2);       
line([size(piSpkMtx,1)/2 size(piSpkMtx,1)/2], get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
line([size(poSpkMtx,1)/2+size(piSpkMtx,1) size(poSpkMtx,1)/2+size(piSpkMtx,1)], get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
line([size(piSpkMtx,1)/2+size(piSpkMtx,1)+size(poSpkMtx,1) size(piSpkMtx,1)/2+size(piSpkMtx,1)+size(poSpkMtx,1)], get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
line([size(poSpkMtx,1)/2+size(piSpkMtx,1)*2+size(poSpkMtx,1) size(poSpkMtx,1)/2+size(piSpkMtx,1)*2+size(poSpkMtx,1)], get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);  
line([size(piSpkMtx,1)/2+size(piSpkMtx,1)*2+size(poSpkMtx,1)*2 size(piSpkMtx,1)/2+size(piSpkMtx,1)*2+size(poSpkMtx,1)*2], get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
line([size(poSpkMtx,1)/2+size(piSpkMtx,1)*3+size(poSpkMtx,1)*2 size(poSpkMtx,1)/2+size(piSpkMtx,1)*3+size(poSpkMtx,1)*2], get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);       
line([size(piSpkMtx,1)/2+size(piSpkMtx,1)*3+size(poSpkMtx,1)*3 size(piSpkMtx,1)/2+size(piSpkMtx,1)*3+size(poSpkMtx,1)*3], get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
line([size(poSpkMtx,1)/2+size(piSpkMtx,1)*4+size(poSpkMtx,1)*3 size(poSpkMtx,1)/2+size(piSpkMtx,1)*4+size(poSpkMtx,1)*3], get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);     
set(gca, 'xtick', []);
ylabel(gca,'Lag (s)');

annotation('textbox', 'position', [0 0.935 1 0.05], 'String', ['\bf\fontsize{10}' sprintf('Bin = %i ms; Step = %i ms; PokeInWindow = [%ims +%ims], PokeOutWindow = [%ims +%ims]', binSize, dsRate, piWindow(1)*1000,piWindow(2)*1000,poWindow(1)*1000, poWindow(2)*1000)],...
    'linestyle', 'none', 'horizontalalignment', 'right');
curDir = cd;
annotation('textbox', 'position', [0.025 0.025 0.7 0.05], 'String', curDir,...
    'linestyle', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');

%% Investigate decoding more.... possibly do rank order quantification for each index... tabulate at each index which odor is 1st, 2nd, 3rd, 4th
% Possibly the indexing will demonstrate a lag effect, i.e. 1st = current odor, 2nd = next odor.
% Alternatively, examine decoding of peaks?
% Alternatively (and I like this one the most), calculate total posterior for different odors, i.e. calculate sum or average posterior value for each odor at each index

%% Step through and decode each OutSeq Trial
outSeqTrials = find(trialInfo(:,6)~=0 & trialInfo(:,5)==1 & trialInfo(:,3)<=4 & trialInfo(:,4)<=4);

osPosts = nan(size(piSpkMtx,1)+size(poSpkMtx,1),(size(piSpkMtx,1)+size(poSpkMtx,1))*2,length(outSeqTrials));
for osT = 1:length(outSeqTrials)
    tempObsv = [piSpkMtx(:,:,outSeqTrials(osT)); poSpkMtx(:,:,outSeqTrials(osT))];
    osPos = trialInfo(outSeqTrials(osT),3);
    osOdr = trialInfo(outSeqTrials(osT),4);
    
    tempPosLikes = [mean(piSpkMtx(:,:,tempISSs(osPos,:)),3); mean(poSpkMtx(:,:,tempISSs(osPos,:)),3)];
    tempOdrLikes = [mean(piSpkMtx(:,:,tempISSs(osOdr,:)),3); mean(poSpkMtx(:,:,tempISSs(osOdr,:)),3)];
    tempLikes = [tempPosLikes; tempOdrLikes];
    
    osPosts(:,:,osT) = CalcStaticBayesPost(tempLikes, tempObsv, binSize);
end
odorLog = [ones(1,size(piSpkMtx,1)+size(poSpkMtx,1))*1,...
    ones(1,size(piSpkMtx,1)+size(poSpkMtx,1))*2];
timeLog = [piTrialTime', poTrialTime'+1+dsRate/1000,...
    piTrialTime', poTrialTime'+1+dsRate/1000];

figure;
subplot(6,4,1:16)
% imagesc(nanmedian(osPosts,3), [0 0.01]); set(gca, 'ydir', 'normal'); colormap jet
imagesc(nanmean(osPosts,3)', [0 0.01]); set(gca, 'ydir', 'normal'); colormap jet
xlabel('Observed');
ylabel('Decodeded');
set(gca,...
    'xtick', [size(piSpkMtx,1)/2, size(poSpkMtx,1)/2+size(piSpkMtx,1)], 'xticklabel', [{'PokeIn'}, {'PokeOut'}],...
    'ytick', [size(piSpkMtx,1)/2 size(piSpkMtx,1) size(piSpkMtx,1)+size(poSpkMtx,1)/2 size(piSpkMtx,1)/2+size(piSpkMtx,1)+size(poSpkMtx,1) size(piSpkMtx,1)*2+size(poSpkMtx,1) size(poSpkMtx,1)/2+size(piSpkMtx,1)*2+size(poSpkMtx,1)],...
    'yticklabel', [{'PokeIn'}, {'\bf\fontsize{14}Position'}, {'PokeOut'}, {'PokeIn'}, {'\bf\fontsize{14}Odor'}, {'PokeOut'}]);
line(get(gca, 'xlim'), [size(piSpkMtx,1) + size(poSpkMtx,1), size(piSpkMtx,1) + size(poSpkMtx,1)], 'linestyle', '-', 'color', 'white', 'linewidth', 2); 
line([size(piSpkMtx,1)/2 size(piSpkMtx,1)/2], get(gca, 'ylim'), 'linestyle', ':', 'color', 'white', 'linewidth', 1);
line([size(poSpkMtx,1)/2+size(piSpkMtx,1) size(poSpkMtx,1)/2+size(piSpkMtx,1)], get(gca, 'ylim'), 'linestyle', ':', 'color', 'white', 'linewidth', 1);

line(get(gca, 'xlim'), [size(piSpkMtx,1)/2 size(piSpkMtx,1)/2], 'linestyle', ':', 'color', 'white', 'linewidth', 1);
line(get(gca, 'xlim'), [size(poSpkMtx,1)/2+size(piSpkMtx,1) size(poSpkMtx,1)/2+size(piSpkMtx,1)], 'linestyle', ':', 'color', 'white', 'linewidth', 1);
line(get(gca, 'xlim'), [size(piSpkMtx,1)/2+size(piSpkMtx,1)*2 size(piSpkMtx,1)/2+size(piSpkMtx,1)*2], 'linestyle', ':', 'color', 'white', 'linewidth', 1);
line(get(gca, 'xlim'), [size(piSpkMtx,1)/2+size(piSpkMtx,1)*2+size(poSpkMtx,1) size(piSpkMtx,1)/2+size(piSpkMtx,1)*2+size(poSpkMtx,1)], 'linestyle', ':', 'color', 'white', 'linewidth', 1);

% Odor Decoding
decodeOdor = DecodeBayesPost(osPosts, odorLog);
posDecode = mean(decodeOdor==1,2);
odrDecode = mean(decodeOdor==2,2);
subplot(6,4,17:20)
plot(1:size(decodeOdor,1), posDecode-odrDecode, 'color', 'k', 'linewidth', 1);
axis tight
set(gca, 'ylim', [-1 1], 'xtick', [size(piSpkMtx,1)/2, size(poSpkMtx,1)/2+size(piSpkMtx,1)], 'xticklabel', [{'PokeIn'}, {'PokeOut'}]);
ylabel([{'Decoded Difference'};{'(Pos-Odor)'}]);
line(get(gca, 'xlim'), [0 0], 'color', 'k', 'linestyle', '--');
line([size(piSpkMtx,1) size(piSpkMtx,1)], get(gca, 'ylim'), 'linestyle', '-', 'color', 'k', 'linewidth', 2);
line([size(piSpkMtx,1)/2 size(piSpkMtx,1)/2], get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
line([size(poSpkMtx,1)/2+size(piSpkMtx,1) size(poSpkMtx,1)/2+size(piSpkMtx,1)], get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);

% Temporal Decoding
decodeTime = DecodeBayesPost(osPosts, timeLog);
decodeLag = nan(size(decodeTime));
for c = 1:size(decodeTime,2)
    decodeLag(:,c) = decodeTime(:,c)-[piTrialTime', poTrialTime'+1+dsRate/1000]';
end
lagMean = nanmean(decodeLag,2);
lagVar = nanstd(decodeLag,1,2);
subplot(6,4,21:24)
plot(lagMean, '-k', 'linewidth', 1);
patch('YData', [lagMean+lagVar; flipud(lagMean-lagVar)],...
    'XData', [1:length(lagMean), length(lagMean):-1:1], 'FaceColor', 'k', 'FaceAlpha', .3, 'edgecolor', 'none');
axis tight
set(gca, 'xtick', [size(piSpkMtx,1)/2, size(poSpkMtx,1)/2+size(piSpkMtx,1)], 'xticklabel', [{'PokeIn'}, {'PokeOut'}]);
line(get(gca, 'xlim'), [0 0], 'color', 'k', 'linestyle', '--');
line([size(piSpkMtx,1) size(piSpkMtx,1)], get(gca, 'ylim'), 'linestyle', '-', 'color', 'k', 'linewidth', 2);
line([size(piSpkMtx,1)/2 size(piSpkMtx,1)/2], get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
line([size(poSpkMtx,1)/2+size(piSpkMtx,1) size(poSpkMtx,1)/2+size(piSpkMtx,1)], get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);

annotation('textbox', 'position', [0 0.935 1 0.05], 'String', ['\bf\fontsize{10}' sprintf('Bin = %i ms; Step = %i ms; PokeInWindow = [%ims +%ims], PokeOutWindow = [%ims +%ims]', binSize, dsRate, piWindow(1)*1000,piWindow(2)*1000,poWindow(1)*1000, poWindow(2)*1000)],...
    'linestyle', 'none', 'horizontalalignment', 'right');
curDir = cd;
annotation('textbox', 'position', [0.025 0.025 0.7 0.05], 'String', curDir,...
    'linestyle', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');

%% Trial After OutSeq
% For TAO, use FIS as likelihoods and TAO as observations
% How to compare TAO with FIS decoding... Compare both posteriors as well
% as the likelihood of decoding the correct odor.
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
        binnedEnsembleMtx(:,u,t) = conv(ensembleMtx(:,u,t), ones(1,binSize)./(binSize/1000), 'same');
%         binnedEnsembleMtx(:,u,t) = conv(trialEnsembleMtx(:,u,t), ones(1,binSize), 'same');
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