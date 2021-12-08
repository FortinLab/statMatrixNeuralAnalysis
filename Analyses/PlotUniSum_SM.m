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
    figure;
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
    set(wire1, 'ylim', [-0.2 0.2], 'ytick', -0.2:0.2:0.2, 'color', 'none', 'xticklabels', []);
    wire1.XAxis.Color = 'none';
    
    wire2 = axes('Position', [0.1 0.85 0.0325 0.05]);
    plot(wire2,xVals,uniSum.TemplateMean{2}, 'linewidth', 1.5, 'color', 'black');
    hold on;
    plot(wire2,xVals,uniSum.TemplateMean{2}+uniSum.TemplateStDev{2}, 'linewidth', 1, 'color', 'black', 'linestyle', ':');
    plot(wire2,xVals,uniSum.TemplateMean{2}-uniSum.TemplateStDev{2}, 'linewidth', 1, 'color', 'black', 'linestyle', ':');
    box off
    axis off
    
    wire3 = axes('Position', [0.15 0.85 0.0325 0.05]);
    plot(wire3,xVals,uniSum.TemplateMean{3}, 'linewidth', 1.5, 'color', 'black');
    hold on;
    plot(wire3,xVals,uniSum.TemplateMean{3}+uniSum.TemplateStDev{3}, 'linewidth', 1, 'color', 'black', 'linestyle', ':');
    plot(wire3,xVals,uniSum.TemplateMean{3}-uniSum.TemplateStDev{3}, 'linewidth', 1, 'color', 'black', 'linestyle', ':');
    box off
    axis off
    
    wire4 = axes('Position', [0.2 0.85 0.0325 0.05]);
    plot(wire4,xVals,uniSum.TemplateMean{4}, 'linewidth', 1.5, 'color', 'black');
    hold on;
    plot(wire4,xVals,uniSum.TemplateMean{4}+uniSum.TemplateStDev{4}, 'linewidth', 1, 'color', 'black', 'linestyle', ':');
    plot(wire4,xVals,uniSum.TemplateMean{4}-uniSum.TemplateStDev{4}, 'linewidth', 1, 'color', 'black', 'linestyle', ':');
    set(wire4, 'xlim', [1 32], 'color', 'none')
    box off
    axis off
    
    linkaxes([wire1, wire2, wire3, wire4], 'xy');
    
    set(wire1, 'ylim', [-0.2 0.2]);
    
    maxTemplate = eval(sprintf('wire%i', uniSum.Spike_Features(1,1)));
    text(uniSum.Spike_Features(2,1), uniSum.Spike_Features(2,2),...
        {'\uparrow', '\fontsize{10} Valley'}, 'horizontalalignment', 'center', 'verticalalignment', 'top', 'parent', maxTemplate);
    text(uniSum.Spike_Features(3,1), uniSum.Spike_Features(3,2),...
        {'\fontsize{10} Peak', '\downarrow'}, 'horizontalalignment', 'center', 'verticalalignment', 'bottom', 'parent', maxTemplate);
    
    %% Plot Mean Evoked Responses
    meanFR = axes('Position', [0.05 0.625 0.135 0.15]);
    BarPlotErrorbars([uniSum.TrialEpochFRs.PreTrialFRSANSA(1), uniSum.TrialEpochFRs.EarlyTrialFRSANSA(1), uniSum.TrialEpochFRs.LateTrialFRSANSA(1), uniSum.TrialEpochFRs.PostTrialFRSANSA(1)],...
        [uniSum.TrialEpochFRs.PreTrialFRSANSA(2), uniSum.TrialEpochFRs.EarlyTrialFRSANSA(2), uniSum.TrialEpochFRs.LateTrialFRSANSA(2), uniSum.TrialEpochFRs.PostTrialFRSANSA(2)]);
    axis(meanFR, 'tight');
    set(meanFR, 'xlim', [0 5], 'xtick', 1:4, 'xticklabel', {'Pre-Trial', 'Early-Trial', 'Late-Trial', 'Post-Trial'}, 'XTickLabelRotation', 45, 'color', 'none');
    box off
%     meanFR.XAxis.Color = 'none';
    title('Trial Epoch')
    ylabel('\bf\fontsize{10}Mean Firing Rate');
    
    fdbkFR = axes('Position', [0.2 0.625 0.06 0.15]);
    BarPlotErrorbars([uniSum.TrialEpochFRs.RewardFRSANSA(1), uniSum.TrialEpochFRs.ErrorFRSANSA(1)],...
        [uniSum.TrialEpochFRs.RewardFRSANSA(2), uniSum.TrialEpochFRs.ErrorFRSANSA(2)]);
    axis(fdbkFR, 'tight');
    set(fdbkFR, 'xlim', [0 3], 'xtick', 1:2, 'xticklabel', {'Reward', 'Error'}, 'XTickLabelRotation', 45, 'color', 'none');
    box off
%     fdbkFR.XAxis.Color = 'none';
    title('Feedback')  
    linkaxes([meanFR fdbkFR], 'y');
    
    curYlim = get(meanFR, 'ylim');
    set(meanFR, 'ylim', [curYlim(1) curYlim(2)+2]);
    curYlim = get(meanFR, 'ylim');
        
    if uniSum.TrialEpochStats.TrialEpochsFSANSA(2)<newCritF
        line(meanFR, [1 4], repmat(curYlim(2)-1, [1 2]), 'color', 'k');
        text(meanFR,2.5, curYlim(2)-1+0.1, '\fontsize{20}\bf*');
    end
    if uniSum.TrialEpochStats.TrialPeriodsFSANSA(2)<newCritF
        line(meanFR, [2 3], repmat(curYlim(2)-2, [1 2]), 'color', 'k');
        text(meanFR, 2.5, curYlim(2)-2+0.1, '\fontsize{20}\bf*');
    end
    if uniSum.TrialEpochStats.FeedbackFSANSA(2)<newCritF
        line(fdbkFR, [1 2], repmat(curYlim(2)-1, [1 2]), 'color', 'k');
        text(fdbkFR, 1.5, curYlim(2)-1+0.1, '\fontsize{20}\bf*');
    end
        

    
    %% Plot Correlation between epochs
    rPcrit = 0.05/14;
    curEpochCorrTable = uniSum.TrialEpochStats.EpochCorrelations.R;
    curEpochCorrTable(5:6,:) = [];
    curEpochSigTable = uniSum.TrialEpochStats.EpochCorrelations.P;
    curEpochSigTable(5:6,:) = [];
    epochCorr = axes('position', [0.065 0.35 0.2 0.175]);
    [r,c] = ind2sub(size(curEpochCorrTable), find(isnan(curEpochCorrTable)));
    imagesc(1:6, 1:4,curEpochCorrTable, [-1 1]);
    hold on;
    for p = 1:length(r)
        patch([c(p)-0.5 c(p)-0.5 c(p)+0.5 c(p)+0.5], [r(p)-0.5 r(p)+0.5 r(p)+0.5 r(p)-0.5], 'white', 'edgecolor', 'white');
    end
    set(gca, 'xtick', 1:6, 'xticklabels', {'\fontsize{8}Pre-Trial', '\fontsize{8}Early-Trial', '\fontsize{8}Late-Trial', '\fontsize{8}Post-Trial', '\fontsize{8}Reward', '\fontsize{8}Error'}, 'XTickLabelRotation', 45,...
        'ytick', 1:4, 'yticklabels', {'\fontsize{8}Pre-Trial', '\fontsize{8}Early-Trial', '\fontsize{8}Late-Trial', '\fontsize{8}Post-Trial'});
    hold on;
    [row,col] = ind2sub(size(curEpochCorrTable),find(curEpochSigTable<rPcrit));
    for roW = 1:length(row)
        text(epochCorr, col(roW), row(roW), '*', 'horizontalalignment', 'center', 'fontweight', 'bold', 'fontsize', 25, 'color', 'k')
    end
    colormap jet
    title('Epoch Correlations');
    cb = colorbar;
    title(cb, '\bfr-Val')
%     drawnow
    
    %% Plot Trial Epoch F-Values
    preTrl = axes('position', [0.05 0.1 0.035 0.15]);    
    bar(preTrl, 1:4, [uniSum.TrialEpochStats.PreTrial.Odor(1),...
        uniSum.TrialEpochStats.PreTrial.OdorSANSA(1),...
        uniSum.TrialEpochStats.PreTrial.Position(1),...
        uniSum.TrialEpochStats.PreTrial.PositionSANSA(1)]);    
    box(preTrl, 'off');
    set(preTrl, 'color', 'none', 'xlim', [0 5], 'xtick', 0.5:3.5,...
        'xticklabel', {'\fontsize{5}Odor', '\fontsize{5}Odor(-A)', '\fontsize{5}Pos', '\fontsize{5}Pos(-A)'}, 'xticklabelrotation', 45);
    ylabel('\bf\fontsize{10}F-Ratio');
    title(preTrl, 'Pre')
    if uniSum.TrialEpochStats.PreTrial.Odor(2)<newCritF
        text(preTrl,1, uniSum.TrialEpochStats.PreTrial.Odor(1)+0.1, '\fontsize{15}\bf*', 'horizontalalignment', 'center');
    end
    if uniSum.TrialEpochStats.PreTrial.OdorSANSA(2)<newCritF
        text(preTrl,2, uniSum.TrialEpochStats.PreTrial.OdorSANSA(1)+0.1, '\fontsize{15}\bf*', 'horizontalalignment', 'center');
    end
    if uniSum.TrialEpochStats.PreTrial.Position(2)<newCritF
        text(preTrl,3, uniSum.TrialEpochStats.PreTrial.Position(1)+0.1, '\fontsize{15}\bf*', 'horizontalalignment', 'center');
    end
    if uniSum.TrialEpochStats.PreTrial.PositionSANSA(2)<newCritF
        text(preTrl,4, uniSum.TrialEpochStats.PreTrial.PositionSANSA(1)+0.1, '\fontsize{15}\bf*', 'horizontalalignment', 'center');
    end
    
    rlyTrl = axes('position', [0.09 0.1 0.035 0.15]);    
    bar(rlyTrl, 1:4, [uniSum.TrialEpochStats.EarlyTrial.Odor(1),...
        uniSum.TrialEpochStats.EarlyTrial.OdorSANSA(1),...
        uniSum.TrialEpochStats.EarlyTrial.Position(1),...
        uniSum.TrialEpochStats.EarlyTrial.PositionSANSA(1)]);
    box(rlyTrl, 'off');    
    set(rlyTrl, 'color', 'none', 'xlim', [0 5], 'yticklabels', [], 'xtick', 0.5:3.5,...
        'xticklabel', {'\fontsize{5}Odor', '\fontsize{5}Odor(-A)', '\fontsize{5}Pos', '\fontsize{5}Pos(-A)'}, 'xticklabelrotation', 45);
    title(rlyTrl, 'Early')
    if uniSum.TrialEpochStats.EarlyTrial.Odor(2)<newCritF
        text(rlyTrl,1, uniSum.TrialEpochStats.EarlyTrial.Odor(1)+1, '\fontsize{15}\bf*', 'horizontalalignment', 'center');
    end
    if uniSum.TrialEpochStats.EarlyTrial.OdorSANSA(2)<newCritF
        text(rlyTrl,2, uniSum.TrialEpochStats.EarlyTrial.OdorSANSA(1)+1, '\fontsize{15}\bf*', 'horizontalalignment', 'center');
    end
    if uniSum.TrialEpochStats.EarlyTrial.Position(2)<newCritF
        text(rlyTrl,3, uniSum.TrialEpochStats.EarlyTrial.Position(1)+1, '\fontsize{15}\bf*', 'horizontalalignment', 'center');
    end
    if uniSum.TrialEpochStats.EarlyTrial.PositionSANSA(2)<newCritF
        text(rlyTrl,4, uniSum.TrialEpochStats.EarlyTrial.PositionSANSA(1)+1, '\fontsize{15}\bf*', 'horizontalalignment', 'center');
    end
    
    ltTrl = axes('position', [0.13 0.1 0.035 0.15]);    
    bar(ltTrl, 1:4, [uniSum.TrialEpochStats.LateTrial.Odor(1),...
        uniSum.TrialEpochStats.LateTrial.OdorSANSA(1),...
        uniSum.TrialEpochStats.LateTrial.Position(1),...
        uniSum.TrialEpochStats.LateTrial.PositionSANSA(1)]);
    box(ltTrl, 'off');
    set(ltTrl, 'color', 'none', 'xlim', [0 5], 'yticklabels', [], 'xtick', 0.5:3.5,...
        'xticklabel', {'\fontsize{5}Odor', '\fontsize{5}Odor(-A)', '\fontsize{5}Pos', '\fontsize{5}Pos(-A)'}, 'xticklabelrotation', 45);
    title(ltTrl, 'Late')     
    if uniSum.TrialEpochStats.LateTrial.Odor(2)<newCritF
        text(ltTrl,1, uniSum.TrialEpochStats.LateTrial.Odor(1)+1, '\fontsize{15}\bf*', 'horizontalalignment', 'center');
    end
    if uniSum.TrialEpochStats.LateTrial.OdorSANSA(2)<newCritF
        text(ltTrl,2, uniSum.TrialEpochStats.LateTrial.OdorSANSA(1)+1, '\fontsize{15}\bf*', 'horizontalalignment', 'center');
    end
    if uniSum.TrialEpochStats.LateTrial.Position(2)<newCritF
        text(ltTrl,3, uniSum.TrialEpochStats.LateTrial.Position(1)+1, '\fontsize{15}\bf*', 'horizontalalignment', 'center');
    end
    if uniSum.TrialEpochStats.LateTrial.PositionSANSA(2)<newCritF
        text(ltTrl,4, uniSum.TrialEpochStats.LateTrial.PositionSANSA(1)+1, '\fontsize{15}\bf*', 'horizontalalignment', 'center');
    end   
    
    pstTrl = axes('position', [0.17 0.1 0.035 0.15]);    
    bar(pstTrl, 1:4, [uniSum.TrialEpochStats.PostTrial.Odor(1),...
        uniSum.TrialEpochStats.PostTrial.OdorSANSA(1),...
        uniSum.TrialEpochStats.PostTrial.Position(1),...
        uniSum.TrialEpochStats.PostTrial.PositionSANSA(1)]);
    box(pstTrl, 'off');
    set(pstTrl, 'color', 'none', 'xlim', [0 5], 'yticklabels', [], 'xtick', 0.5:3.5,...
        'xticklabel', {'\fontsize{5}Odor', '\fontsize{5}Odor(-A)', '\fontsize{5}Pos', '\fontsize{5}Pos(-A)'}, 'xticklabelrotation', 45);    
    title(pstTrl, 'Post')         
    if uniSum.TrialEpochStats.PostTrial.Odor(2)<newCritF
        text(pstTrl,1, uniSum.TrialEpochStats.PostTrial.Odor(1)+1, '\fontsize{15}\bf*', 'horizontalalignment', 'center');
    end
    if uniSum.TrialEpochStats.PostTrial.OdorSANSA(2)<newCritF
        text(pstTrl,2, uniSum.TrialEpochStats.PostTrial.OdorSANSA(1)+1, '\fontsize{15}\bf*', 'horizontalalignment', 'center');
    end
    if uniSum.TrialEpochStats.PostTrial.Position(2)<newCritF
        text(pstTrl,3, uniSum.TrialEpochStats.PostTrial.Position(1)+1, '\fontsize{15}\bf*', 'horizontalalignment', 'center');
    end
    if uniSum.TrialEpochStats.PostTrial.PositionSANSA(2)<newCritF
        text(pstTrl,4, uniSum.TrialEpochStats.PostTrial.PositionSANSA(1)+1, '\fontsize{15}\bf*', 'horizontalalignment', 'center');
    end       
    
    rwd = axes('position', [0.21 0.1 0.035 0.15]);    
    bar(rwd, 1:4, [uniSum.TrialEpochStats.Reward.Odor(1),...
        uniSum.TrialEpochStats.Reward.OdorSANSA(1),...
        uniSum.TrialEpochStats.Reward.Position(1),...
        uniSum.TrialEpochStats.Reward.PositionSANSA(1)]);
    box(rwd, 'off');
    set(rwd, 'color', 'none', 'xlim', [0 5], 'yticklabels', [], 'xtick', 0.5:3.5,...
        'xticklabel', {'\fontsize{5}Odor', '\fontsize{5}Odor(-A)', '\fontsize{5}Pos', '\fontsize{5}Pos(-A)'}, 'xticklabelrotation', 45);            
    title(rwd, 'Rwd')               
    if uniSum.TrialEpochStats.Reward.Odor(2)<newCritF
        text(rwd,1, uniSum.TrialEpochStats.Reward.Odor(1)+1, '\fontsize{15}\bf*', 'horizontalalignment', 'center');
    end
    if uniSum.TrialEpochStats.Reward.OdorSANSA(2)<newCritF
        text(rwd,2, uniSum.TrialEpochStats.Reward.OdorSANSA(1)+1, '\fontsize{15}\bf*', 'horizontalalignment', 'center');
    end
    if uniSum.TrialEpochStats.Reward.Position(2)<newCritF
        text(rwd,3, uniSum.TrialEpochStats.Reward.Position(1)+1, '\fontsize{15}\bf*', 'horizontalalignment', 'center');
    end
    if uniSum.TrialEpochStats.Reward.PositionSANSA(2)<newCritF
        text(rwd,4, uniSum.TrialEpochStats.Reward.PositionSANSA(1)+1, '\fontsize{15}\bf*', 'horizontalalignment', 'center');
    end       
    
    err = axes('position', [0.25 0.1 0.035 0.15]);    
    bar(err, 1:4, [uniSum.TrialEpochStats.Error.Odor(1),...
        uniSum.TrialEpochStats.Error.OdorSANSA(1),...
        uniSum.TrialEpochStats.Error.Position(1),...
        uniSum.TrialEpochStats.Error.PositionSANSA(1)]);
    box(err, 'off'); 
    set(err, 'color', 'none', 'xlim', [0 5], 'yticklabels', [], 'xtick', 0.5:1.5,...
        'xticklabel', {'\fontsize{5}Odor', '\fontsize{5}Odor(-A)', '\fontsize{5}Pos', '\fontsize{5}Pos(-A)'}, 'xticklabelrotation', 45);            
    title(err, 'Err')                           
    if uniSum.TrialEpochStats.Error.Odor(2)<newCritF
        text(err,1, uniSum.TrialEpochStats.Error.Odor(1)+1, '\fontsize{15}\bf*', 'horizontalalignment', 'center');
    end
%     if uniSum.TrialEpochStats.Error.OdorSANSA(2)<newCritF
%         text(err,2, uniSum.TrialEpochStats.Error.OdorSANSA(1)+1, '\fontsize{15}\bf*', 'horizontalalignment', 'center');
%     end
    if uniSum.TrialEpochStats.Error.Position(2)<newCritF
        text(err,3, uniSum.TrialEpochStats.Error.Position(1)+1, '\fontsize{15}\bf*', 'horizontalalignment', 'center');
    end
    if uniSum.TrialEpochStats.Error.PositionSANSA(2)<newCritF
        text(err,4, uniSum.TrialEpochStats.Error.PositionSANSA(1)+1, '\fontsize{15}\bf*', 'horizontalalignment', 'center');
    end                   
    
    linkaxes([preTrl, rlyTrl, ltTrl, pstTrl, rwd, err], 'y');
    
    %% Plot Trial Rasters Aligned Poke In
    corIStrlLog = uniSum.TrialInfo.InSeqTrialLog & uniSum.TrialInfo.CorrectTrialLog;
%     corIStrlLog = uniSum.TrialInfo.CorrectTrialLog;
    pokeInTrialTimeBins = uniSum.WholeTrial.PokeIn.TimeBins;
    pokeDur = uniSum.TrialInfo.PokeDuration;
    pokeDur(~corIStrlLog) = [];   
    positionVals = uniSum.TrialInfo.Position;
    positionVals(~corIStrlLog) = [];
    rasters = uniSum.WholeTrial.PokeIn.Rasters;
    rasters(~corIStrlLog) = []; 
    
    posSortingMtx = [positionVals; 1:length(positionVals); pokeDur]';
    sortedPosMtx = sortrows(posSortingMtx);
    
    rasters = rasters(sortedPosMtx(:,2))';
    trlPokeDur = nan(size(rasters));
    trlPokeDurTrlNum = nan(size(rasters));
    rasterTrl = cell(size(rasters));
    rasterTrlPos = cell(size(rasters));
    for r = 1:length(rasters)
        rasterTrl{r} = repmat(r, size(rasters{r}));
        rasterTrlPos{r} = repmat(sortedPosMtx(r,1), size(rasters{r}));
        if ~isempty(rasters{r})
            trlPokeDur(r) = sortedPosMtx(r,3);
            trlPokeDurTrlNum(r) = r;
        end
    end
    pokeInRasterScatVals = [cell2mat(rasterTrlPos), cell2mat(rasterTrl), cell2mat(rasters)];
    trlPokeDur(isnan(trlPokeDur)) = [];
    trlPokeDurTrlNum(isnan(trlPokeDurTrlNum)) = [];
    
    pokeInScatPlot = axes('Position', [0.325 0.7 0.3 0.25]);
    hold on;
    aLog = pokeInRasterScatVals(:,1)==1;
    scatter(pokeInRasterScatVals(aLog,3), pokeInRasterScatVals(aLog,2)*-1, 10, [44/255 168/255 224/255], 'filled');
    bLog = pokeInRasterScatVals(:,1)==2;
    scatter(pokeInRasterScatVals(bLog,3), pokeInRasterScatVals(bLog,2)*-1, 10, [154/255 133/255 122/255], 'filled');
    cLog = pokeInRasterScatVals(:,1)==3;
    scatter(pokeInRasterScatVals(cLog,3), pokeInRasterScatVals(cLog,2)*-1, 10, [9/255 161/255 74/255], 'filled');
    dLog = pokeInRasterScatVals(:,1)==4;
    scatter(pokeInRasterScatVals(dLog,3), pokeInRasterScatVals(dLog,2)*-1, 10, [128/255 66/255 151/255], 'filled');
    axis(pokeInScatPlot, 'tight');
    pokeInScatPlot.YAxis.Color = 'none';
    scatter(trlPokeDur, trlPokeDurTrlNum*-1, 20, 'k', '<', 'filled');
    set(pokeInScatPlot, 'xlim', [min(pokeInTrialTimeBins) 0.5], 'color', 'none');
    line([0 0], get(pokeInScatPlot, 'ylim'), 'color','k');
    title('PokeIn Aligned');
    
    % Plot Poke Out Aligned
        pokeOutTrialTimeBins = uniSum.WholeTrial.PokeOut.TimeBins;
    pokeDur = uniSum.TrialInfo.PokeDuration;
    pokeDur(~corIStrlLog) = [];
    positionVals = uniSum.TrialInfo.Position;
    positionVals(~corIStrlLog) = [];
    pokeOutRasters = uniSum.WholeTrial.PokeOut.Rasters;
    pokeOutRasters(~corIStrlLog) = [];

    posSortingMtx = [positionVals; 1:length(positionVals); pokeDur]';
    sortedPosMtx = sortrows(posSortingMtx);

    pokeOutRasters = pokeOutRasters(sortedPosMtx(:,2))';
    trlPokeDur = nan(size(pokeOutRasters));
    trlPokeDurTrlNum = nan(size(pokeOutRasters));
    rasterTrl = cell(size(pokeOutRasters));
    rasterTrlPos = cell(size(pokeOutRasters));
    for r = 1:length(pokeOutRasters)
        rasterTrl{r} = repmat(r, size(pokeOutRasters{r}));
        rasterTrlPos{r} = repmat(sortedPosMtx(r,1), size(pokeOutRasters{r}));
        if ~isempty(pokeOutRasters{r})
            trlPokeDur(r) = sortedPosMtx(r,3);
            trlPokeDurTrlNum(r) = r;
        end
    end
    pokeOutSasterScatVals = [cell2mat(rasterTrlPos), cell2mat(rasterTrl), cell2mat(pokeOutRasters)];
    trlPokeDur(isnan(trlPokeDur)) = [];
    trlPokeDurTrlNum(isnan(trlPokeDurTrlNum)) = [];

    pokeOutScatPlot = axes('Position', [0.65 0.7 0.3 0.25]);
    hold on;
    aLog = pokeOutSasterScatVals(:,1)==1;
    scatter(pokeOutSasterScatVals(aLog,3), pokeOutSasterScatVals(aLog,2)*-1, 10, [44/255 168/255 224/255], 'filled');
    bLog = pokeOutSasterScatVals(:,1)==2;
    scatter(pokeOutSasterScatVals(bLog,3), pokeOutSasterScatVals(bLog,2)*-1, 10, [154/255 133/255 122/255], 'filled');
    cLog = pokeOutSasterScatVals(:,1)==3;
    scatter(pokeOutSasterScatVals(cLog,3), pokeOutSasterScatVals(cLog,2)*-1, 10, [9/255 161/255 74/255], 'filled');
    dLog = pokeOutSasterScatVals(:,1)==4;
    scatter(pokeOutSasterScatVals(dLog,3), pokeOutSasterScatVals(dLog,2)*-1, 10, [128/255 66/255 151/255], 'filled');
    axis(pokeOutScatPlot, 'tight');
    pokeOutScatPlot.YAxis.Color = 'none';
    scatter(trlPokeDur*-1, trlPokeDurTrlNum*-1, 20, 'k', '>', 'filled');
    set(pokeOutScatPlot, 'xlim', [-0.5 max(pokeOutTrialTimeBins)], 'color', 'none');
    line([0 0], get(pokeOutScatPlot, 'ylim'), 'color','k');
    title('PokeOut Aligned');
    
    
    %% Plot Trial Evoked Activity       
    corIStrlLog = uniSum.TrialInfo.InSeqTrialLog & uniSum.TrialInfo.CorrectTrialLog;
    aLog = uniSum.TrialInfo.Position==1 & corIStrlLog;
    bLog = uniSum.TrialInfo.Position==2 & corIStrlLog;
    cLog = uniSum.TrialInfo.Position==3 & corIStrlLog;
    dLog = uniSum.TrialInfo.Position==4 & corIStrlLog;
            
    pokeInTrialTimeBins = uniSum.WholeTrial.PokeIn.TimeBins;
    pokeInFRplot = axes('Position', [0.325 0.4 0.3 0.25]);
    hold on;
    pokeInFR = uniSum.WholeTrial.PokeIn.FiringRate;
    pokeInA = pokeInFR(:,aLog);
    pokeInB = pokeInFR(:,bLog);
    pokeInC = pokeInFR(:,cLog);
    pokeInD = pokeInFR(:,dLog);
    patch(pokeInFRplot, 'XData', [pokeInTrialTimeBins', fliplr(pokeInTrialTimeBins')],...
        'YData', [(mean(pokeInA,2)+(std(pokeInA,1,2)./sqrt(size(pokeInA,2)-1)))', fliplr((mean(pokeInA,2)-(std(pokeInA,1,2)./sqrt(size(pokeInA,2)-1)))')],...
        'FaceColor', 'none', 'FaceAlpha', 0.15,...
        'EdgeColor', [44/255 168/255 224/255]);
    patch(pokeInFRplot, 'XData', [pokeInTrialTimeBins', fliplr(pokeInTrialTimeBins')],...
        'YData', [(mean(pokeInB,2)+(std(pokeInB,1,2)./sqrt(size(pokeInB,2)-1)))', fliplr((mean(pokeInB,2)-(std(pokeInB,1,2)./sqrt(size(pokeInB,2)-1)))')],...
        'FaceColor', 'none', 'FaceAlpha', 0.15,...
        'EdgeColor', [154/255 133/255 122/255]); 
    patch(pokeInFRplot, 'XData', [pokeInTrialTimeBins', fliplr(pokeInTrialTimeBins')],...
        'YData', [(mean(pokeInC,2)+(std(pokeInC,1,2)./sqrt(size(pokeInC,2)-1)))', fliplr((mean(pokeInC,2)-(std(pokeInC,1,2)./sqrt(size(pokeInC,2)-1)))')],...
        'FaceColor', 'none','FaceAlpha', 0.15,...
        'EdgeColor', [9/255 161/255 74/255]);
    patch(pokeInFRplot, 'XData', [pokeInTrialTimeBins', fliplr(pokeInTrialTimeBins')],...
        'YData', [(mean(pokeInD,2)+(std(pokeInD,1,2)./sqrt(size(pokeInD,2)-1)))', fliplr((mean(pokeInD,2)-(std(pokeInD,1,2)./sqrt(size(pokeInD,2)-1)))')],...
        'FaceColor', 'none', 'FaceAlpha', 0.15,...
        'EdgeColor', [128/255 66/255 151/255]);
    hold on;
    plot(pokeInFRplot,pokeInTrialTimeBins, mean(pokeInA,2), 'color', [44/255 168/255 224/255], 'linewidth', 2);
    plot(pokeInFRplot,pokeInTrialTimeBins, mean(pokeInB,2), 'color', [154/255 133/255 122/255], 'linewidth', 2);    
    plot(pokeInFRplot,pokeInTrialTimeBins, mean(pokeInC,2), 'color', [9/255 161/255 74/255], 'linewidth', 2);
    plot(pokeInFRplot,pokeInTrialTimeBins, mean(pokeInD,2), 'color', [128/255 66/255 151/255], 'linewidth', 2);
    set(pokeInFRplot, 'xlim', [min(pokeInTrialTimeBins) 0.5], 'color', 'none');
    line([0 0], get(pokeInFRplot, 'ylim'), 'color','k');
    ylabel('\bfFiring Rate (spk/s)');
        
    pokeOutTrialTimeBins = uniSum.WholeTrial.PokeOut.TimeBins;
    pokeOutFRplot = axes('Position', [0.65 0.4 0.3 0.25]);
    hold on;
    pokeOutFR = uniSum.WholeTrial.PokeOut.FiringRate;
    pokeOutA = pokeOutFR(:,aLog);
    pokeOutB = pokeOutFR(:,bLog);   
    pokeOutC = pokeOutFR(:,cLog);
    pokeOutD = pokeOutFR(:,dLog);
    
    hold on;
    patch(pokeOutFRplot, 'XData', [pokeOutTrialTimeBins', fliplr(pokeOutTrialTimeBins')],...
        'YData', [(mean(pokeOutA,2)+(std(pokeOutA,1,2)./sqrt(size(pokeOutA,2)-1)))', fliplr((mean(pokeOutA,2)-(std(pokeOutA,1,2)./sqrt(size(pokeOutA,2)-1)))')],...
        'FaceColor', 'none', 'FaceAlpha', 0.15,...
        'EdgeColor', [44/255 168/255 224/255]); 
    patch(pokeOutFRplot, 'XData', [pokeOutTrialTimeBins', fliplr(pokeOutTrialTimeBins')],...
        'YData', [(mean(pokeOutB,2)+(std(pokeOutB,1,2)./sqrt(size(pokeOutB,2)-1)))', fliplr((mean(pokeOutB,2)-(std(pokeOutB,1,2)./sqrt(size(pokeOutB,2)-1)))')],...
        'FaceColor', 'none', 'FaceAlpha', 0.15,...
        'EdgeColor', [154/255 133/255 122/255]); 
    patch(pokeOutFRplot, 'XData', [pokeOutTrialTimeBins', fliplr(pokeOutTrialTimeBins')],...
        'YData', [(mean(pokeOutC,2)+(std(pokeOutC,1,2)./sqrt(size(pokeOutC,2)-1)))', fliplr((mean(pokeOutC,2)-(std(pokeOutC,1,2)./sqrt(size(pokeOutC,2)-1)))')],...
        'FaceColor', 'none','FaceAlpha', 0.15,...
        'EdgeColor', [9/255 161/255 74/255]);
    patch(pokeOutFRplot, 'XData', [pokeOutTrialTimeBins', fliplr(pokeOutTrialTimeBins')],...
        'YData', [(mean(pokeOutD,2)+(std(pokeOutD,1,2)./sqrt(size(pokeOutD,2)-1)))', fliplr((mean(pokeOutD,2)-(std(pokeOutD,1,2)./sqrt(size(pokeOutD,2)-1)))')],...
        'FaceColor', 'none', 'FaceAlpha', 0.15,...
        'EdgeColor', [128/255 66/255 151/255]);
    aPlot = plot(pokeOutFRplot,pokeOutTrialTimeBins, mean(pokeOutA,2), 'color', [44/255 168/255 224/255], 'linewidth', 2);
    bPlot = plot(pokeOutFRplot,pokeOutTrialTimeBins, mean(pokeOutB,2), 'color', [154/255 133/255 122/255], 'linewidth', 2);
    cPlot = plot(pokeOutFRplot,pokeOutTrialTimeBins, mean(pokeOutC,2), 'color', [9/255 161/255 74/255], 'linewidth', 2);
    dPlot = plot(pokeOutFRplot,pokeOutTrialTimeBins, mean(pokeOutD,2), 'color', [128/255 66/255 151/255], 'linewidth', 2);
    set(pokeOutFRplot, 'xlim', [-0.5 max(pokeOutTrialTimeBins)], 'color', 'none');
    line([0 0], get(pokeOutFRplot, 'ylim'), 'color','k');
    pokeOutFRplot.YAxis.Color = 'none';
    leg = legend([aPlot, bPlot, cPlot, dPlot], {'A', 'B', 'C', 'D'});
    leg.Orientation = 'horizontal';
    legendPos = leg.Position;
    leg.Position = [0.58 0.65 legendPos(3) legendPos(4)]; 
    
    linkaxes([pokeInFRplot, pokeOutFRplot], 'y');
    yLimVals = get(pokeInFRplot, 'ylim');
    set(pokeInFRplot, 'ylim', [0 max(yLimVals)]);
        
    %% Plot F-Ratios for Position and Odor
    pokeInTrialTimeBins = uniSum.WholeTrial.PokeIn.TimeBins;
    pokeInPosFZval = uniSum.InformationContent.TrialPokeIn.PosZ;
    pokeInPosSANSAfZval = uniSum.InformationContent.TrialPokeInSANSA.PosZ;
    pokeInOdrFZval = uniSum.InformationContent.TrialPokeIn.OdorZ;
    pokeInOdrSANSAfZval = uniSum.InformationContent.TrialPokeInSANSA.OdorZ;
    pokeInPrevOdrSANSAfZval = uniSum.InformationContent.TrialPokeInSANSA.PrevOdorZ;

    pokeInFZplot = axes('Position', [0.325 0.1 0.3 0.25]);
    plot(pokeInTrialTimeBins, pokeInPosFZval, 'linewidth', 1, 'color', 'k');
    hold on;
    plot(pokeInTrialTimeBins, pokeInOdrFZval, 'linewidth', 1, 'color', 'r');
    plot(pokeInTrialTimeBins, pokeInPosSANSAfZval, 'linewidth', 1, 'color', 'k', 'linestyle', ':');
%     hold on;
    plot(pokeInTrialTimeBins, pokeInOdrSANSAfZval, 'linewidth', 1, 'color', 'r', 'linestyle', ':');
    plot(pokeInTrialTimeBins, pokeInPrevOdrSANSAfZval, 'linewidth', 1, 'color', 'c', 'linestyle', ':');
    line(pokeInFZplot,[min(pokeInTrialTimeBins) max(pokeInTrialTimeBins)], [0 0], 'linewidth', 1.5, 'color', 'k');
    line(pokeInFZplot,[min(pokeInTrialTimeBins) max(pokeInTrialTimeBins)], [2 2], 'linewidth', 0.5, 'color', 'k', 'linestyle', '-.');
    line(pokeInFZplot,[min(pokeInTrialTimeBins) max(pokeInTrialTimeBins)], [-2 -2], 'linewidth', 0.5, 'color', 'k', 'linestyle', '-.');
    set(pokeInFZplot, 'xlim', [min(pokeInTrialTimeBins) 0.5], 'color', 'none');
    %     pokeInFZplot.XAxis.Color = 'none';
    ylabel({'\bfF-Ratios(z-norm)'});
    box off

    pokeOutTrialTimeBins = uniSum.WholeTrial.PokeOut.TimeBins;
    pokeOutPosFZval = uniSum.InformationContent.TrialPokeOut.PosZ;
    pokeOutPosSANSAfZval = uniSum.InformationContent.TrialPokeOutSANSA.PosZ;
    pokeOutOdrFZval = uniSum.InformationContent.TrialPokeOut.OdorZ;
    pokeOutOdrSANSAfZval = uniSum.InformationContent.TrialPokeOutSANSA.OdorZ;
    pokeOutPrevOdrSANSAfZval = uniSum.InformationContent.TrialPokeOutSANSA.PrevOdorZ;

    pokeOutFZplot = axes('Position', [0.65 0.1 0.3 0.25]);
    posPlot = plot(pokeOutTrialTimeBins, pokeOutPosFZval, 'linewidth', 1, 'color', 'k');
    hold on;
    odrPlot = plot(pokeOutTrialTimeBins, pokeOutOdrFZval, 'linewidth', 1, 'color', 'r');
    posSANSAplot = plot(pokeOutTrialTimeBins, pokeOutPosSANSAfZval, 'linewidth', 1, 'color', 'k', 'linestyle', ':');
    hold on
    odrSANSAplot = plot(pokeOutTrialTimeBins, pokeOutOdrSANSAfZval, 'linewidth', 1, 'color', 'r', 'linestyle', ':');
    prevOdrSANSAplot = plot(pokeInTrialTimeBins, pokeOutPrevOdrSANSAfZval, 'linewidth', 1, 'color', 'c', 'linestyle', ':');
    line(pokeOutFZplot,[min(pokeOutTrialTimeBins) max(pokeOutTrialTimeBins)], [0 0], 'linewidth', 1.5, 'color', 'k');
    line(pokeOutFZplot,[min(pokeOutTrialTimeBins) max(pokeOutTrialTimeBins)], [2 2], 'linewidth', 0.5, 'color', 'k', 'linestyle', '-.');
    line(pokeOutFZplot,[min(pokeOutTrialTimeBins) max(pokeOutTrialTimeBins)], [-2 -2], 'linewidth', 0.5, 'color', 'k', 'linestyle', '-.');
    set(pokeOutFZplot, 'xlim', [-0.5 max(pokeOutTrialTimeBins)], 'color', 'none');
%         pokeOutFZplot.XAxis.Color = 'none';
    box off
    leg = legend([posPlot, odrPlot, posSANSAplot, odrSANSAplot, prevOdrSANSAplot], {'Position', 'Odor', 'Pos(sansA)', 'Odor(sansA)', 'PrevOdor'});
%     leg = legend([posSANSAplot, odrSANSAplot, prevOdrSANSAplot], 'Pos(sansA)', 'Odor(sansA)', 'PrevOdor');
    leg.Orientation = 'horizontal';
    legendPos = leg.Position;
    leg.Position = [0.4 0.35 legendPos(3) legendPos(4)];
    
    linkaxes([pokeInFZplot, pokeOutFZplot], 'y');
    line(pokeInFZplot,[0 0], get(pokeInFZplot, 'ylim'), 'color','k');
    line(pokeOutFZplot,[0 0], get(pokeOutFZplot, 'ylim'), 'color','k');
    
    %% Print/Save output
    drawnow
    orient(gcf, 'tall');
    orient(gcf, 'landscape');
%     print
    print('-painters', gcf, '-dpdf', uniSumFiles{u}(1:end-4));
    saveas(gcf, [uniSumFiles{u}(1:end-4) '.png'], 'png');
    close gcf;
end
    