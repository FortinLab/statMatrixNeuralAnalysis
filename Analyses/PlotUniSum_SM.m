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
    
    %% 
    uniSum.TrialEpochStats.PreTrial
    

    