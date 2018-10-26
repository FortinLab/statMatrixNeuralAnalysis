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
%% Run Analysis
for u = 1:length(uniSumFiles)
    load(uniSumFiles{u});
    
    %% Plot Waveform Template
    wire1 = axes('Position', [0.05 0.9 0.0325 0.05]);
    xVals = 1:length(uniSum.TemplateMean{1});
    plot(wire1,xVals,uniSum.TemplateMean{1}, 'linewidth', 1.5, 'color', 'black');
    hold on;
    plot(wire1,xVals,uniSum.TemplateMean{1}+uniSum.TemplateStDev{1}, 'linewidth', 1, 'color', 'black', 'linestyle', ':');
    plot(wire1,xVals,uniSum.TemplateMean{1}-uniSum.TemplateStDev{1}, 'linewidth', 1, 'color', 'black', 'linestyle', ':');
    box off
%     title({[uniSum.UnitName]; 'Wire1'});
%     title('Wire1')
    
    wire2 = axes('Position', [0.1 0.9 0.0325 0.05]);
    plot(wire2,xVals,uniSum.TemplateMean{2}, 'linewidth', 1.5, 'color', 'black');
    hold on;
    plot(wire2,xVals,uniSum.TemplateMean{2}+uniSum.TemplateStDev{2}, 'linewidth', 1, 'color', 'black', 'linestyle', ':');
    plot(wire2,xVals,uniSum.TemplateMean{2}-uniSum.TemplateStDev{2}, 'linewidth', 1, 'color', 'black', 'linestyle', ':');
    box off
%     title({[uniSum.UnitName]; 'Wire2'});
%     title('Wire2')
    
    wire3 = axes('Position', [0.15 0.9 0.0325 0.05]);
    plot(wire3,xVals,uniSum.TemplateMean{3}, 'linewidth', 1.5, 'color', 'black');
    hold on;
    plot(wire3,xVals,uniSum.TemplateMean{3}+uniSum.TemplateStDev{3}, 'linewidth', 1, 'color', 'black', 'linestyle', ':');
    plot(wire3,xVals,uniSum.TemplateMean{3}-uniSum.TemplateStDev{3}, 'linewidth', 1, 'color', 'black', 'linestyle', ':');
    box off
%     title({[uniSum.UnitName]; 'Wire3'});
%     title('Wire3')
    
    wire4 = axes('Position', [0.2 0.9 0.0325 0.05]);
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
    meanFR = axes('Position', [0.05 0.68 0.135 0.15]);
    BarPlotErrorbars([uniSum.TrialEpochFRs.PreTrialFR(1), uniSum.TrialEpochFRs.EarlyTrialFR(1), uniSum.TrialEpochFRs.LateTrialFR(1), uniSum.TrialEpochFRs.PostTrialFR(1)],...
        [uniSum.TrialEpochFRs.PreTrialFR(2), uniSum.TrialEpochFRs.EarlyTrialFR(2), uniSum.TrialEpochFRs.LateTrialFR(2), uniSum.TrialEpochFRs.PostTrialFR(2)]);
    axis(meanFR, 'tight');
    set(meanFR, 'xlim', [0 5], 'xtick', 1:4, 'xticklabel', {'Pre-Trial', 'Early-Trial', 'Late-Trial', 'Post-Trial'}, 'XTickLabelRotation', 45, 'color', 'none');
    box off
    title('Mean Epoch Activity')
    
    fdbkFR = axes('Position', [0.2 0.68 0.06 0.15]);
    BarPlotErrorbars([uniSum.TrialEpochFRs.RewardFR(1), uniSum.TrialEpochFRs.ErrorFR(1)],...
        [uniSum.TrialEpochFRs.RewardFR(2), uniSum.TrialEpochFRs.ErrorFR(2)]);
    axis(fdbkFR, 'tight');
    set(fdbkFR, 'xlim', [0 3], 'xtick', 1:2, 'xticklabel', {'Reward', 'Error'}, 'XTickLabelRotation', 45, 'color', 'none');
    box off
    title('Feedback Activity')
    
    
    linkaxes([meanFR fdbkFR], 'y');
    
    
     
     
     
     
    
    
%% Summary 1: 