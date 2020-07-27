smPath = uigetdir;
cd(smPath);
files = dir(smPath);
fileNames = {files.name};
% Identify list of all statMatrix files
decodingFiles = fileNames(cellfun(@(a)~isempty(a), regexp(fileNames, 'Decodings')))';

%%
decodes = cell(1,length(decodingFiles),4);
osDiffs = cell(1,length(decodingFiles));
for f = 1:length(decodingFiles)
    load(decodingFiles{f});
    decodes{1,f,1} = decodings.Odor(:,1);
    decodes{1,f,2} = decodings.Odor(:,2);
    decodes{1,f,3} = decodings.Odor(:,3);
    decodes{1,f,4} = decodings.Odor(:,4);
    osDiffs{f} = decodings.OSdiff;
end
decodes = cell2mat(decodes);
osDiffs = cell2mat(osDiffs);

figure;
plot(1:size(decodes,1), nanmean(decodes(:,:,1),2), 'color', [44/255 168/255 224/255], 'linewidth', 1);
hold on;
plot(1:size(decodes,1), nanmean(decodes(:,:,2),2), 'color', [154/255 133/255 122/255], 'linewidth', 1);
plot(1:size(decodes,1), nanmean(decodes(:,:,3),2), 'color', [9/255 161/255 74/255], 'linewidth', 1);
plot(1:size(decodes,1), nanmean(decodes(:,:,4),2), 'color', [128/255 66/255 151/255], 'linewidth', 1);
legend('A', 'B', 'C','D', 'location', 'southoutside', 'orientation', 'horizontal');
for o = 1:4
    tempPostMean = nanmean(decodes(:,:,o),2);
%     tempPostVar = nanstd(decodes(:,:,o),1,2);
    tempPostVar = SEMcalc(decodes(:,:,o)')';
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
ylabel([{'Decoded Odor'};{'% Trials'}]);
xlabel('Sequence Time');
title('GE11 Decodings (All InSeq, Leave-1-Out; Mean +/- SEM; 5xSessions)');

%%
figure; 
plot(1:size(osDiffs,1), nanmean(osDiffs,2), 'color', 'k', 'linewidth', 1);
axis tight
patch('YData', [nanmean(osDiffs,2)+SEMcalc(osDiffs')'; flipud(nanmean(osDiffs,2)-SEMcalc(osDiffs')')],...
    'XData', [1:length(nanmean(osDiffs,2)), length(nanmean(osDiffs,2)):-1:1], 'FaceAlpha', .3);
set(gca, 'ylim', [-1 1], 'xtick', [size(decodings.PItime,1)/2, size(decodings.POtime,1)/2+size(decodings.PItime,1)], 'xticklabel', [{'PokeIn'}, {'PokeOut'}]);
ylabel([{'Decoded Difference'};{'(Pos-Odor)'}]);
line(get(gca, 'xlim'), [0 0], 'color', 'k', 'linestyle', '--');
line([size(decodings.PItime,1) size(decodings.PItime,1)], get(gca, 'ylim'), 'linestyle', '-', 'color', 'k', 'linewidth', 2);
line([size(decodings.PItime,1)/2 size(decodings.PItime,1)/2], get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
line([size(decodings.POtime,1)/2+size(decodings.POtime,1) size(decodings.POtime,1)/2+size(decodings.PItime,1)], get(gca, 'ylim'), 'linestyle', ':', 'color', 'k', 'linewidth', 2);
title('GE11 OutSeq Decoding Diff (Position-Odor; Mean +/- SEM; 5xSessions)');


