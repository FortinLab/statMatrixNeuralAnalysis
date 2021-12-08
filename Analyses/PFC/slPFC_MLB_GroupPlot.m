smPath = uigetdir;
cd(smPath);
files = dir(smPath);
fileNames = {files.name};
% Identify list of all statMatrix files
decodingFiles = fileNames(cellfun(@(a)~isempty(a), regexp(fileNames, 'Decodings')))';

%%
decodes = cell(1,length(decodingFiles),4);
osDiffs = cell(1,length(decodingFiles));
skpDiffs = cell(1,length(decodingFiles));
repDiffs = cell(1,length(decodingFiles));
trlPrdDecodes = cell(1,length(decodingFiles));
ssnSMI = nan(1,length(decodingFiles));
ssnSMIsfp = nan(1,length(decodingFiles));
trlPrdDprms = nan(4,length(decodingFiles));
for f = 1:length(decodingFiles)
    load(decodingFiles{f});
    decodes{1,f,1} = decodings.Odor(:,1);
    decodes{1,f,2} = decodings.Odor(:,2);
    decodes{1,f,3} = decodings.Odor(:,3);
    decodes{1,f,4} = decodings.Odor(:,4);
    osDiffs{f} = decodings.OSdiff;
    skpDiffs{f} = decodings.OSdiffSKP;
    repDiffs{f} = decodings.OSdiffREP;
    trlPrdDecodes{f} = decodings.TrialPeriods;
    ssnSMI(f) = decodings.SMIbehav;
    ssnSMIsfp(f) = decodings.SMIbehavSFP;
    trlPrdDprms(:,f) = decodings.TrialPeriodsDprm;
end
decodes = cell2mat(decodes);
osDiffs = cell2mat(osDiffs);
skpDiffs = cell2mat(skpDiffs);
repDiffs = cell2mat(repDiffs);

%%
figure;
plot(1:size(decodes,1), smooth(nanmean(decodes(:,:,1),2)), 'color', [44/255 168/255 224/255], 'linewidth', 1);
hold on;
plot(1:size(decodes,1), smooth(nanmean(decodes(:,:,2),2)), 'color', [154/255 133/255 122/255], 'linewidth', 1);
plot(1:size(decodes,1), smooth(nanmean(decodes(:,:,3),2)), 'color', [9/255 161/255 74/255], 'linewidth', 1);
plot(1:size(decodes,1), smooth(nanmean(decodes(:,:,4),2)), 'color', [128/255 66/255 151/255], 'linewidth', 1);
legend('1', '2', '3','4', 'location', 'southoutside', 'orientation', 'horizontal');
for o = 1:4
    tempPostMean = smooth(nanmean(decodes(:,:,o),2));
%     tempPostVar = nanstd(decodes(:,:,o),1,2);
    tempPostVar = smooth(SEMcalc(decodes(:,:,o)')');
    tp = patch('YData', [tempPostMean+tempPostVar; flipud(tempPostMean-tempPostVar)],...
        'XData', [1:length(tempPostMean), length(tempPostMean):-1:1], 'FaceAlpha', .3);
    if o==1
        set(tp, 'FaceColor', [44/255 168/255 224/255], 'edgecolor', [44/255 168/255 224/255], 'EdgeAlpha', .5);
    elseif o==2
        set(tp, 'FaceColor', [154/255 133/255 122/255], 'edgecolor', [154/255 133/255 122/255], 'EdgeAlpha', .5);
    elseif o==3
        set(tp, 'FaceColor', [9/255 161/255 74/255], 'edgecolor', [9/255 161/255 74/255], 'EdgeAlpha', .5);
    else
        set(tp, 'FaceColor', [128/255 66/255 151/255], 'edgecolor', [128/255 66/255 151/255], 'EdgeAlpha', .5);
    end 
end

axis tight
set(gca, 'ylim', [0 1], 'xtick', []);
line([size(decodings.PItime,1) + size(decodings.POtime,1), size(decodings.PItime,1) + size(decodings.POtime,1)], get(gca, 'ylim'), 'linestyle', '-', 'color', 'k', 'linewidth', 2);   
line([size(decodings.PItime,1) + size(decodings.POtime,1), size(decodings.PItime,1) + size(decodings.POtime,1)]*2, get(gca, 'ylim'), 'linestyle', '-', 'color', 'k', 'linewidth', 2);       
line([size(decodings.PItime,1) + size(decodings.POtime,1), size(decodings.PItime,1) + size(decodings.POtime,1)]*3, get(gca, 'ylim'), 'linestyle', '-', 'color', 'k', 'linewidth', 2);       
line([size(decodings.PItime,1) + size(decodings.POtime,1), size(decodings.PItime,1) + size(decodings.POtime,1)]*4, get(gca, 'ylim'), 'linestyle', '-', 'color', 'k', 'linewidth', 2);       
line([size(decodings.PItime,1)/2 size(decodings.PItime,1)/2], get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
line([size(decodings.POtime,1)/2+size(decodings.PItime,1) size(decodings.POtime,1)/2+size(decodings.PItime,1)], get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
line([size(decodings.PItime,1)/2+size(decodings.PItime,1)+size(decodings.POtime,1) size(decodings.PItime,1)/2+size(decodings.PItime,1)+size(decodings.POtime,1)], get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
line([size(decodings.POtime,1)/2+size(decodings.PItime,1)*2+size(decodings.POtime,1) size(decodings.POtime,1)/2+size(decodings.PItime,1)*2+size(decodings.POtime,1)], get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);  
line([size(decodings.PItime,1)/2+size(decodings.PItime,1)*2+size(decodings.POtime,1)*2 size(decodings.PItime,1)/2+size(decodings.PItime,1)*2+size(decodings.POtime,1)*2], get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
line([size(decodings.POtime,1)/2+size(decodings.PItime,1)*3+size(decodings.POtime,1)*2 size(decodings.POtime,1)/2+size(decodings.PItime,1)*3+size(decodings.POtime,1)*2], get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);       
line([size(decodings.PItime,1)/2+size(decodings.PItime,1)*3+size(decodings.POtime,1)*3 size(decodings.PItime,1)/2+size(decodings.PItime,1)*3+size(decodings.POtime,1)*3], get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
line([size(decodings.POtime,1)/2+size(decodings.PItime,1)*4+size(decodings.POtime,1)*3 size(decodings.POtime,1)/2+size(decodings.PItime,1)*4+size(decodings.POtime,1)*3], get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);     
ylabel([{'Decoded Position'};{'% Trials'}]);
xlabel('Sequence Time');
title(sprintf('GE11 Decodings (All InSeq, Leave-1-Out; Mean +/- SEM; %ixSessions)', length(decodingFiles)));

%%
odorLog = [ones(1,size(decodings.PItime,1)+size(decodings.POtime,1))*1,...
    ones(1,size(decodings.PItime,1)+size(decodings.POtime,1))*2,...
    ones(1,size(decodings.PItime,1)+size(decodings.POtime,1))*3,...
    ones(1,size(decodings.PItime,1)+size(decodings.POtime,1))*4];

figure;
sps = nan(1,4);
for p = 1:4
    decLog = odorLog==p;
    sps(p) = subplot(1,4,p);
    plot(1:sum(decLog), smooth(nanmean(decodes(decLog,:,1),2)), 'color', [44/255 168/255 224/255], 'linewidth', 1);
    hold on;
    plot(1:sum(decLog), smooth(nanmean(decodes(decLog,:,2),2)), 'color', [154/255 133/255 122/255], 'linewidth', 1);
    plot(1:sum(decLog), smooth(nanmean(decodes(decLog,:,3),2)), 'color', [9/255 161/255 74/255], 'linewidth', 1);
    plot(1:sum(decLog), smooth(nanmean(decodes(decLog,:,4),2)), 'color', [128/255 66/255 151/255], 'linewidth', 1);
    for o = 1:4
        tempPostMean = smooth(nanmean(decodes(decLog,:,o),2));
        tempPostVar = smooth(SEMcalc(decodes(decLog,:,o)')');
        tp = patch('YData', [tempPostMean+tempPostVar; flipud(tempPostMean-tempPostVar)],...
            'XData', [1:length(tempPostMean), length(tempPostMean):-1:1], 'FaceAlpha', .3);
        if o==1
            set(tp, 'FaceColor', [44/255 168/255 224/255], 'edgecolor', [44/255 168/255 224/255], 'EdgeAlpha', .5);
        elseif o==2
            set(tp, 'FaceColor', [154/255 133/255 122/255], 'edgecolor', [154/255 133/255 122/255], 'EdgeAlpha', .5);
        elseif o==3
            set(tp, 'FaceColor', [9/255 161/255 74/255], 'edgecolor', [9/255 161/255 74/255], 'EdgeAlpha', .5);
        else
            set(tp, 'FaceColor', [128/255 66/255 151/255], 'edgecolor', [128/255 66/255 151/255], 'EdgeAlpha', .5);
        end
    end
    axis(sps(p), 'tight');
    set(sps(p), 'ylim', [0 1], 'xticklabels', []);
    line([size(decodings.PItime,1)/2 size(decodings.PItime,1)/2], get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
    line([size(decodings.POtime,1)/2+size(decodings.PItime,1) size(decodings.POtime,1)/2+size(decodings.PItime,1)], get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
end
    legend('1', '2', '3','4', 'location', 'southoutside', 'orientation', 'horizontal');

%%
figure; 
plot(1:size(osDiffs,1), smooth(nanmean(osDiffs,2)), 'color', 'k', 'linewidth', 1);
axis tight
patch('YData', [smooth(nanmean(osDiffs,2)+SEMcalc(osDiffs')'); flipud(smooth(nanmean(osDiffs,2)-SEMcalc(osDiffs')'))],...
    'XData', [1:length(nanmean(osDiffs,2)), length(nanmean(osDiffs,2)):-1:1], 'FaceAlpha', .3);
set(gca, 'ylim', [-1 1], 'xtick', [size(decodings.PItime,1)/2, size(decodings.POtime,1)/2+size(decodings.PItime,1)], 'xticklabel', [{'PokeIn'}, {'PokeOut'}]);
ylabel([{'Decoded Difference'};{'(Pos-Odor)'}]);
line(get(gca, 'xlim'), [0 0], 'color', 'k', 'linestyle', '--');
line([size(decodings.PItime,1) size(decodings.PItime,1)], get(gca, 'ylim'), 'linestyle', '-', 'color', 'k', 'linewidth', 2);
line([size(decodings.PItime,1)/2 size(decodings.PItime,1)/2], get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
line([size(decodings.POtime,1)/2+size(decodings.PItime,1) size(decodings.POtime,1)/2+size(decodings.PItime,1)], get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
title(sprintf('GE11 OutSeq Decoding Diff (Position-Odor; Mean +/- SEM; %ixSessions)', length(decodingFiles)));

%%
figure;
skpPlot = plot(1:size(skpDiffs,1), smooth(nanmean(skpDiffs,2)), 'color', 'k', 'linewidth', 1);
patch('YData', [smooth(nanmean(skpDiffs,2)+SEMcalc(skpDiffs')'); smooth(flipud(nanmean(skpDiffs,2)-SEMcalc(skpDiffs')'))],...
    'XData', [1:length(nanmean(skpDiffs,2)), length(nanmean(skpDiffs,2)):-1:1], 'FaceAlpha', .3);
axis tight
hold on;
repPlot = plot(1:size(repDiffs,1), smooth(nanmean(repDiffs,2)), 'color', 'r', 'linewidth', 1);
patch('YData', [smooth(nanmean(repDiffs,2)+SEMcalc(repDiffs')'); smooth(flipud(nanmean(repDiffs,2)-SEMcalc(repDiffs')'))],...
    'XData', [1:length(nanmean(repDiffs,2)), length(nanmean(repDiffs,2)):-1:1], 'FaceColor', 'r', 'EdgeColor', 'r', 'FaceAlpha', .3);
set(gca, 'ylim', [-1 1], 'xtick', [size(decodings.PItime,1)/2, size(decodings.POtime,1)/2+size(decodings.PItime,1)], 'xticklabel', [{'PokeIn'}, {'PokeOut'}]);
legend([skpPlot repPlot], 'Skips', 'Repeats');
ylabel([{'Decoded Difference'};{'(Pos-Odor)'}]);
line(get(gca, 'xlim'), [0 0], 'color', 'k', 'linestyle', '--');
line([size(decodings.PItime,1) size(decodings.PItime,1)], get(gca, 'ylim'), 'linestyle', '-', 'color', 'k', 'linewidth', 2);
line([size(decodings.PItime,1)/2 size(decodings.PItime,1)/2], get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
line([size(decodings.POtime,1)/2+size(decodings.PItime,1) size(decodings.POtime,1)/2+size(decodings.PItime,1)], get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
title('OutSeq Diff by Type');

%%
meanTrialDecodes = nan(4,4,4);
for p = 1:4
    meanTrialDecodes(:,:,p) = mean(cell2mat(reshape(cellfun(@(a)a(:,:,p),trlPrdDecodes, 'uniformoutput',0), [1,1,length(trlPrdDecodes)])),3);
end

figure;
subplot(2,4,1)
imagesc(meanTrialDecodes(:,:,1), [0 0.5]);
set(gca, 'xtick', 1:4, 'ytick', 1:4);
title('Pre-Trial Period');
ylabel('Decoded Position');
subplot(2,4,2)
imagesc(meanTrialDecodes(:,:,2), [0 0.5]);
set(gca, 'xtick', 1:4, 'ytick', 1:4);
title('Early-Trial Period');
subplot(2,4,3)
imagesc(meanTrialDecodes(:,:,3), [0 0.5]);
set(gca, 'xtick', 1:4, 'ytick', 1:4);
title('Late-Trial Period');
xlabel('True Position');
ylabel('Decoded Position');
subplot(2,4,4)
imagesc(meanTrialDecodes(:,:,4), [0 0.5]);
set(gca, 'xtick', 1:4, 'ytick', 1:4);
title('Post-Trial Period');
xlabel('True Position');
colormap jet

subplot(2,4,5)
corrScatPlot(ssnSMI', trlPrdDprms(1,:)', 'SMI', 'dPrm');
subplot(2,4,6)
corrScatPlot(ssnSMI', trlPrdDprms(2,:)', 'SMI', 'dPrm');
subplot(2,4,7)
corrScatPlot(ssnSMI', trlPrdDprms(3,:)', 'SMI', 'dPrm');
subplot(2,4,8)
corrScatPlot(ssnSMI', trlPrdDprms(4,:)', 'SMI', 'dPrm');