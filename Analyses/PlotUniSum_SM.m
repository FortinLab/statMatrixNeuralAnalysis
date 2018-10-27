function PlotUniSum_SM
%% PlotUniSum_SM
%   Function to plot the data stored in the unitSummary (suffix '_UniSum')
%   files.

%%
files = dir(cd);
fileNames = {files.name};
matFiles = fileNames(cell2mat(cellfun(@(a)~isempty(a), strfind(fileNames, '.mat'), 'uniformoutput', 0)))';
if isempty(matFiles)
    matDir = uigetdir(cd, 'Select the folder with the statMatrix Files');
    if matDir==0
        disp('Analysis Cancelled')
        return
    else
        cd(matDir)
    end
end
uniSumLog = false(size(matFiles));
for fl = 1:length(matFiles)
    variableInfo = who('-file', matFiles{fl});
    if sum(ismember(variableInfo, 'uniSum'))==1
        uniSumLog(fl) = true;
    end
end
if sum(uniSumLog) == 0
    error('No Unit Summary files located, check directory');
end

uniSumFiles = matFiles(uniSumLog);
%% Declare runtime variables
newCritF = 0.05/23;
%% Run Analysis
for u = 1:length(uniSumFiles)
    load(uniSumFiles{u});
    figure
    %% Title
    dirParts = strsplit(uniSum.Directory, '\');
    figTitle = annotation('textbox', 'position', [0.025 0.935 0.7 0.05], 'String', ['\bf\fontsize{20}' sprintf('%s %s',uniSum.UnitName, dirParts{end})],...
        'linestyle', 'none', 'horizontalalignment', 'left');
    %% Plot Waveform Template
    wire1 = axes('Position', [0.05 0.85 0.0325 0.05]);
    xVals = 1:length(uniSum.TemplateMean{1});
    plot(wire1,xVals,uniSum.TemplateMean{1}, 'linewidth', 1.5, 'color', 'black');
    hold on;
    plot(wire1,xVals,uniSum.TemplateMean{1}+uniSum.TemplateStDev{1}, 'linewidth', 1, 'color', 'black', 'linestyle', ':');
    plot(wire1,xVals,uniSum.TemplateMean{1}-uniSum.TemplateStDev{1}, 'linewidth', 1, 'color', 'black', 'linestyle', ':');
    box off
%     title({[uniSum.UnitName]; 'Wire1'});
%     title('Wire1')
    
    wire2 = axes('Position', [0.125 0.85 0.0325 0.05]);
    plot(wire2,xVals,uniSum.TemplateMean{2}, 'linewidth', 1.5, 'color', 'black');
    hold on;
    plot(wire2,xVals,uniSum.TemplateMean{2}+uniSum.TemplateStDev{2}, 'linewidth', 1, 'color', 'black', 'linestyle', ':');
    plot(wire2,xVals,uniSum.TemplateMean{2}-uniSum.TemplateStDev{2}, 'linewidth', 1, 'color', 'black', 'linestyle', ':');
    box off
%     title({[uniSum.UnitName]; 'Wire2'});
%     title('Wire2')
    
    wire3 = axes('Position', [0.175 0.85 0.0325 0.05]);
    plot(wire3,xVals,uniSum.TemplateMean{3}, 'linewidth', 1.5, 'color', 'black');
    hold on;
    plot(wire3,xVals,uniSum.TemplateMean{3}+uniSum.TemplateStDev{3}, 'linewidth', 1, 'color', 'black', 'linestyle', ':');
    plot(wire3,xVals,uniSum.TemplateMean{3}-uniSum.TemplateStDev{3}, 'linewidth', 1, 'color', 'black', 'linestyle', ':');
    box off
%     title({[uniSum.UnitName]; 'Wire3'});
%     title('Wire3')
    
    wire4 = axes('Position', [0.225 0.85 0.0325 0.05]);
    plot(wire4,xVals,uniSum.TemplateMean{4}, 'linewidth', 1.5, 'color', 'black');
    hold on;
    plot(wire4,xVals,uniSum.TemplateMean{4}+uniSum.TemplateStDev{4}, 'linewidth', 1, 'color', 'black', 'linestyle', ':');
    plot(wire4,xVals,uniSum.TemplateMean{4}-uniSum.TemplateStDev{4}, 'linewidth', 1, 'color', 'black', 'linestyle', ':');
    set(wire4, 'xlim', [1 32], 'color', 'none')
    box off
%     title({[uniSum.UnitName]; 'Wire4'});
%     title('Wire4')
    
    linkaxes([wire1, wire2, wire3, wire4], 'xy');
    
    maxTemplate = eval(sprintf('wire%i', uniSum.Spike_Features(1,1)));
    axis(maxTemplate, 'tight')
    set(maxTemplate, 'color', 'none', 'xticklabel', []);
    
    wireHandleNames = {'wire1', 'wire2', 'wire3', 'wire4'};
    nonMaxTemplate = wireHandleNames(~strcmp(sprintf('wire%i', uniSum.Spike_Features(1,1)), wireHandleNames));
    for t = 1:length(nonMaxTemplate)
        axis(eval(nonMaxTemplate{t}), 'off');
%         set(eval(nonMaxTemplate{t}), 'axis', 'off');
    end
    
    text(uniSum.Spike_Features(2,1), uniSum.Spike_Features(2,2),...
        {'\uparrow', '\fontsize{10} Valley'}, 'horizontalalignment', 'center', 'verticalalignment', 'top', 'parent', maxTemplate);
    text(uniSum.Spike_Features(3,1), uniSum.Spike_Features(3,2),...
        {'\fontsize{10} Peak', '\downarrow'}, 'horizontalalignment', 'center', 'verticalalignment', 'bottom', 'parent', maxTemplate);
    
    %% Plot Mean Evoked Responses
    meanFR = axes('Position', [0.05 0.6 0.135 0.15]);
    BarPlotErrorbars([uniSum.TrialEpochFRs.PreTrialFR(1), uniSum.TrialEpochFRs.EarlyTrialFR(1), uniSum.TrialEpochFRs.LateTrialFR(1), uniSum.TrialEpochFRs.PostTrialFR(1)],...
        [uniSum.TrialEpochFRs.PreTrialFR(2), uniSum.TrialEpochFRs.EarlyTrialFR(2), uniSum.TrialEpochFRs.LateTrialFR(2), uniSum.TrialEpochFRs.PostTrialFR(2)]);
    axis(meanFR, 'tight');
    set(meanFR, 'xlim', [0 5], 'xtick', 1:4, 'xticklabel', {'Pre-Trial', 'Early-Trial', 'Late-Trial', 'Post-Trial'}, 'XTickLabelRotation', 45, 'color', 'none');
    box off
    title('Mean Epoch Activity')
    
    fdbkFR = axes('Position', [0.2 0.6 0.06 0.15]);
    BarPlotErrorbars([uniSum.TrialEpochFRs.RewardFR(1), uniSum.TrialEpochFRs.ErrorFR(1)],...
        [uniSum.TrialEpochFRs.RewardFR(2), uniSum.TrialEpochFRs.ErrorFR(2)]);
    axis(fdbkFR, 'tight');
    set(fdbkFR, 'xlim', [0 3], 'xtick', 1:2, 'xticklabel', {'Reward', 'Error'}, 'XTickLabelRotation', 45, 'color', 'none');
    box off
    title('Feedback Activity')  
    linkaxes([meanFR fdbkFR], 'y');
    
    curYlim = get(meanFR, 'ylim');
    set(meanFR, 'ylim', [curYlim(1) curYlim(2)+2]);
    curYlim = get(meanFR, 'ylim');
        
    if uniSum.TrialEpochStats.TrialEpochsF(2)<newCritF
        line(meanFR, [1 4], repmat(curYlim(2)-1, [1 2]), 'color', 'k');
        text(meanFR,2.5, curYlim(2)-1+0.1, '\fontsize{20}\bf*');
    end
    if uniSum.TrialEpochStats.TrialPeriodsF(2)<newCritF
        line(meanFR, [2 3], repmat(curYlim(2)-2, [1 2]), 'color', 'k');
        text(meanFR, 2.5, curYlim(2)-2+0.1, '\fontsize{20}\bf*');
    end
    if uniSum.TrialEpochStats.FeedbackF(2)<newCritF
        line(fdbkFR, [1 2], repmat(curYlim(2)-1, [1 2]), 'color', 'k');
        text(fdbkFR, 1.5, curYlim(2)-1+0.1, '\fontsize{20}\bf*');
    end
        

    
    %% Plot Correlation between epochs
    rPcrit = 0.05/14;
    curEpochCorrTable = uniSum.TrialEpochStats.EpochCorrelations.R;
    curEpochSigTable = uniSum.TrialEpochStats.EpochCorrelations.P;
    epochCorr = axes('position', [0.065 0.3 0.2 0.2]);
    [r,c] = ind2sub(size(curEpochCorrTable), find(isnan(curEpochCorrTable)));
    imagesc(curEpochCorrTable, [-1 1]);
    hold on;
    for p = 1:length(r)
        patch([c(p)-0.5 c(p)-0.5 c(p)+0.5 c(p)+0.5], [r(p)-0.5 r(p)+0.5 r(p)+0.5 r(p)-0.5], 'white', 'edgecolor', 'white');
    end
    set(gca, 'xticklabels', {'Pre-Trial', 'Early-Trial', 'Late-Trial', 'Post-Trial', 'Reward', 'Error'}, 'XTickLabelRotation', 45,...
        'yticklabels', {'Pre-Trial', 'Early-Trial', 'Late-Trial', 'Post-Trial', 'Reward', 'Error'});
    hold on;
    [col,row] = ind2sub([6 6],find(curEpochSigTable<rPcrit));
    for roW = 1:length(row)
        text(row, col, '*', 'horizontalalignment', 'center', 'fontweight', 'bold', 'fontsize', 20, 'color', 'white')
    end
    colormap jet
    title('Trial Epoch Correlations');
%     colorbar
    cb = colorbar;
    cb.Label.String = '\bf r-Val';
    drawnow
    
    %% Plot Trial Epoch F-Values
    preTrl = axes('position', [0.05 0.05 0.027 0.15]);    
    bar(preTrl, 1:4, [uniSum.TrialEpochStats.PreTrial.PrevOdor(1),...
        uniSum.TrialEpochStats.PreTrial.Odor(1),...
        uniSum.TrialEpochStats.PreTrial.Position(1),...
        uniSum.TrialEpochStats.PreTrial.NxtPos(1)]);    
    box(preTrl, 'off');
    set(preTrl, 'color', 'none', 'xlim', [0 5], 'xtick', 1:4,...
        'xticklabel', {'PrvOdr', 'Odor', 'Pos', 'NxtPos'}, 'xticklabelrotation', 45);
    title(preTrl, 'Pre-Trial')
    if uniSum.TrialEpochStats.PreTrial.PrevOdor(2)<newCritF
        text(preTrl,1, uniSum.TrialEpochStats.PreTrial.PrevOdor(1)+0.1, '\fontsize{15}\bf*', 'horizontalalignment', 'center');
    end
    if uniSum.TrialEpochStats.PreTrial.Odor(2)<newCritF
        text(preTrl,2, uniSum.TrialEpochStats.PreTrial.Odor(1)+0.1, '\fontsize{15}\bf*', 'horizontalalignment', 'center');
    end
    if uniSum.TrialEpochStats.PreTrial.Position(2)<newCritF
        text(preTrl,3, uniSum.TrialEpochStats.PreTrial.Position(1)+0.1, '\fontsize{15}\bf*', 'horizontalalignment', 'center');
    end
    if uniSum.TrialEpochStats.PreTrial.NxtPos(2)<newCritF
        text(preTrl,4, uniSum.TrialEpochStats.PreTrial.NxtPos(1)+0.1, '\fontsize{15}\bf*', 'horizontalalignment', 'center');
    end
    
    rlyTrl = axes('position', [0.09 0.05 0.027 0.15]);    
    bar(rlyTrl, 1:4, [uniSum.TrialEpochStats.EarlyTrial.PrevOdor(1),...
        uniSum.TrialEpochStats.EarlyTrial.Odor(1),...
        uniSum.TrialEpochStats.EarlyTrial.Position(1),...
        uniSum.TrialEpochStats.EarlyTrial.NxtPos(1)]);
    box(rlyTrl, 'off');    
    set(rlyTrl, 'color', 'none', 'xlim', [0 5], 'yticklabels', [], 'xtick', 1:4,...
        'xticklabel', {'PrvOdr', 'Odor', 'Pos', 'NxtPos'}, 'xticklabelrotation', 45);
    title(rlyTrl, 'Early-Trial')
    if uniSum.TrialEpochStats.PreTrial.PrevOdor(2)<newCritF
        text(rlyTrl,1, uniSum.TrialEpochStats.PrevOdor.PrevOdor(1)+1, '\fontsize{15}\bf*', 'horizontalalignment', 'center');
    end
    if uniSum.TrialEpochStats.EarlyTrial.Odor(2)<newCritF
        text(rlyTrl,2, uniSum.TrialEpochStats.EarlyTrial.Odor(1)+1, '\fontsize{15}\bf*', 'horizontalalignment', 'center');
    end
    if uniSum.TrialEpochStats.EarlyTrial.Position(2)<newCritF
        text(rlyTrl,3, uniSum.TrialEpochStats.EarlyTrial.Position(1)+1, '\fontsize{15}\bf*', 'horizontalalignment', 'center');
    end
    if uniSum.TrialEpochStats.EarlyTrial.NxtPos(2)<newCritF
        text(rlyTrl,4, uniSum.TrialEpochStats.EarlyTrial.NxtPos(1)+1, '\fontsize{15}\bf*', 'horizontalalignment', 'center');
    end
    
    ltTrl = axes('position', [0.13 0.05 0.027 0.15]);    
    bar(ltTrl, 1:4, [uniSum.TrialEpochStats.LateTrial.PrevOdor(1),...
        uniSum.TrialEpochStats.LateTrial.Odor(1),...
        uniSum.TrialEpochStats.LateTrial.Position(1),...
        uniSum.TrialEpochStats.LateTrial.NxtPos(1)]);
    box(ltTrl, 'off');
    set(ltTrl, 'color', 'none', 'xlim', [0 5], 'yticklabels', [], 'xtick', 1:4,...
        'xticklabel', {'PrvOdr', 'Odor', 'Pos', 'NxtPos'}, 'xticklabelrotation', 45);
    title(ltTrl, 'Late-Trial')     
    if uniSum.TrialEpochStats.LateTrial.PrevOdor(2)<newCritF
        text(ltTrl,1, uniSum.TrialEpochStats.LateTrial.PrevOdor(1)+1, '\fontsize{15}\bf*', 'horizontalalignment', 'center');
    end
    if uniSum.TrialEpochStats.LateTrial.Odor(2)<newCritF
        text(ltTrl,2, uniSum.TrialEpochStats.LateTrial.Odor(1)+1, '\fontsize{15}\bf*', 'horizontalalignment', 'center');
    end
    if uniSum.TrialEpochStats.LateTrial.Position(2)<newCritF
        text(ltTrl,3, uniSum.TrialEpochStats.LateTrial.Position(1)+1, '\fontsize{15}\bf*', 'horizontalalignment', 'center');
    end
    if uniSum.TrialEpochStats.LateTrial.NxtPos(2)<newCritF
        text(ltTrl,4, uniSum.TrialEpochStats.LateTrial.NxtPos(1)+1, '\fontsize{15}\bf*', 'horizontalalignment', 'center');
    end   
    
    pstTrl = axes('position', [0.17 0.05 0.027 0.15]);    
    bar(pstTrl, 1:4, [uniSum.TrialEpochStats.PostTrial.PrevOdor(1),...
        uniSum.TrialEpochStats.PostTrial.Odor(1),...
        uniSum.TrialEpochStats.PostTrial.Position(1),...
        uniSum.TrialEpochStats.PostTrial.NxtPos(1)]);
    box(pstTrl, 'off');
    set(pstTrl, 'color', 'none', 'xlim', [0 5], 'yticklabels', [], 'xtick', 1:4,...
        'xticklabel', {'PrvOdr', 'Odor', 'Pos', 'NxtPos'}, 'xticklabelrotation', 45);    
    title(pstTrl, 'Post-Trial')         
    if uniSum.TrialEpochStats.PostTrial.PrevOdor(2)<newCritF
        text(pstTrl,1, uniSum.TrialEpochStats.PostTrial.PrevOdor(1)+1, '\fontsize{15}\bf*', 'horizontalalignment', 'center');
    end
    if uniSum.TrialEpochStats.PostTrial.Odor(2)<newCritF
        text(pstTrl,2, uniSum.TrialEpochStats.PostTrial.Odor(1)+1, '\fontsize{15}\bf*', 'horizontalalignment', 'center');
    end
    if uniSum.TrialEpochStats.PostTrial.Position(2)<newCritF
        text(pstTrl,3, uniSum.TrialEpochStats.PostTrial.Position(1)+1, '\fontsize{15}\bf*', 'horizontalalignment', 'center');
    end
    if uniSum.TrialEpochStats.PostTrial.NxtPos(2)<newCritF
        text(pstTrl,4, uniSum.TrialEpochStats.PostTrial.NxtPos(1)+1, '\fontsize{15}\bf*', 'horizontalalignment', 'center');
    end       
    
    rwd = axes('position', [0.21 0.05 0.027 0.15]);    
    bar(rwd, 1:2, [uniSum.TrialEpochStats.Reward.Odor(1),...
        uniSum.TrialEpochStats.Reward.Position(1)]);
    box(rwd, 'off');
    set(rwd, 'color', 'none', 'xlim', [0 3], 'yticklabels', [], 'xtick', 1:2,...
        'xticklabel', {'Odor', 'Pos'}, 'xticklabelrotation', 45);            
    title(rwd, 'Reward')               
    if uniSum.TrialEpochStats.Reward.Odor(2)<newCritF
        text(rwd,1, uniSum.TrialEpochStats.Reward.Odor(1)+1, '\fontsize{15}\bf*', 'horizontalalignment', 'center');
    end
    if uniSum.TrialEpochStats.Reward.Position(2)<newCritF
        text(rwd,2, uniSum.TrialEpochStats.Reward.Position(1)+1, '\fontsize{15}\bf*', 'horizontalalignment', 'center');
    end
    
    err = axes('position', [0.25 0.05 0.027 0.15]);    
    bar(err, 1:2, [uniSum.TrialEpochStats.Error.Odor(1),...
        uniSum.TrialEpochStats.Error.Position(1)]);
    box(err, 'off'); 
    set(err, 'color', 'none', 'xlim', [0 3], 'yticklabels', [], 'xtick', 1:2,...
        'xticklabel', {'Odor', 'Pos'}, 'xticklabelrotation', 45);            
    title(err, 'Error')                           
    if uniSum.TrialEpochStats.Error.Odor(2)<newCritF
        text(err,1, uniSum.TrialEpochStats.Error.Odor(1)+1, '\fontsize{15}\bf*', 'horizontalalignment', 'center');
    end
    if uniSum.TrialEpochStats.Error.Position(2)<newCritF
        text(err,2, uniSum.TrialEpochStats.Error.Position(1)+1, '\fontsize{15}\bf*', 'horizontalalignment', 'center');
    end                   
    
    linkaxes([preTrl, rlyTrl, ltTrl, pstTrl, rwd, err], 'y');
    
    %% Plot Trial Rasters
    trialTimeBins = uniSum.WholeTrial.TimeBins;
    positionVals = uniSum.TrialInfo.Position;
    posSortingMtx = [positionVals; 1:length(positionVals)]';
    sortedPosMtx = sortrows(posSortingMtx);
    
    rasters = uniSum.WholeTrial.Rasters(sortedPosMtx(:,2))';
    rasterTrl = cell(size(rasters));
    rasterTrlPos = cell(size(rasters));
    for r = 1:length(rasters)
        rasterTrl{r} = repmat(r, size(rasters{r}));
        rasterTrlPos{r} = repmat(sortedPosMtx(r,1), size(rasters{r}));
    end
    rasterScatVals = [cell2mat(rasterTrlPos), cell2mat(rasterTrl), cell2mat(rasters)];
    
    scatPlot = axes('Position', [0.3 0.7 0.65 0.25]);
    hold on;
    aLog = rasterScatVals(:,1)==1;
    scatter(rasterScatVals(aLog,3), rasterScatVals(aLog,2), 10, [44/255 168/255 224/255], 'filled');
    bLog = rasterScatVals(:,1)==2;
    scatter(rasterScatVals(bLog,3), rasterScatVals(bLog,2), 10, [154/255 133/255 122/255], 'filled');
    cLog = rasterScatVals(:,1)==3;
    scatter(rasterScatVals(cLog,3), rasterScatVals(cLog,2), 10, [9/255 161/255 74/255], 'filled');
    dLog = rasterScatVals(:,1)==3;
    scatter(rasterScatVals(dLog,3), rasterScatVals(dLog,2), 10, [128/255 66/255 151/255], 'filled');
    set(scatPlot, 'xlim', [min(trialTimeBins) max(trialTimeBins)], 'color', 'none');
    line([0 0], get(scatPlot, 'ylim'), 'color','k');
    axis(scatPlot, 'off');
    
    
    
        
        

    