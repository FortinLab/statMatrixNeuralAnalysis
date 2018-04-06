function PlotUnitFeatures_SM(curUnitSummary, curUniSpikeTimes, figID, saveYN)
%% PlotUnitFeatures
%   Create a figure that shows the features of the unit
%
%%
if nargin==2
    figID = figure;
else
    figure(figID);
end
%% Plot Waveform Template
wire1 = subplot(2,4,1);
xVals = 1:length(curUnitSummary.TemplateMean{1});
plot(wire1,xVals,curUnitSummary.TemplateMean{1}, 'linewidth', 1.5, 'color', 'black');
hold on;
plot(wire1,xVals,curUnitSummary.TemplateMean{1}+curUnitSummary.TemplateStDev{1}, 'linewidth', 1, 'color', 'black', 'linestyle', ':');
plot(wire1,xVals,curUnitSummary.TemplateMean{1}-curUnitSummary.TemplateStDev{1}, 'linewidth', 1, 'color', 'black', 'linestyle', ':');
box off
axis off
title([curUnitSummary.UnitName ' Wire1']);

wire2 = subplot(2,4,2);
plot(wire2,xVals,curUnitSummary.TemplateMean{2}, 'linewidth', 1.5, 'color', 'black');
hold on;
plot(wire2,xVals,curUnitSummary.TemplateMean{2}+curUnitSummary.TemplateStDev{2}, 'linewidth', 1, 'color', 'black', 'linestyle', ':');
plot(wire2,xVals,curUnitSummary.TemplateMean{2}-curUnitSummary.TemplateStDev{2}, 'linewidth', 1, 'color', 'black', 'linestyle', ':');
box off
axis off
title([curUnitSummary.UnitName ' Wire2']);

wire3 = subplot(2,4,3);
plot(wire3,xVals,curUnitSummary.TemplateMean{3}, 'linewidth', 1.5, 'color', 'black');
hold on;
plot(wire3,xVals,curUnitSummary.TemplateMean{3}+curUnitSummary.TemplateStDev{3}, 'linewidth', 1, 'color', 'black', 'linestyle', ':');
plot(wire3,xVals,curUnitSummary.TemplateMean{3}-curUnitSummary.TemplateStDev{3}, 'linewidth', 1, 'color', 'black', 'linestyle', ':');
box off
axis off
title([curUnitSummary.UnitName ' Wire3']);

wire4 = subplot(2,4,4);
plot(wire4,xVals,curUnitSummary.TemplateMean{4}, 'linewidth', 1.5, 'color', 'black');
hold on;
plot(wire4,xVals,curUnitSummary.TemplateMean{4}+curUnitSummary.TemplateStDev{4}, 'linewidth', 1, 'color', 'black', 'linestyle', ':');
plot(wire4,xVals,curUnitSummary.TemplateMean{4}-curUnitSummary.TemplateStDev{4}, 'linewidth', 1, 'color', 'black', 'linestyle', ':');
set(wire4, 'xlim', [1 32])
box off
axis off
title([curUnitSummary.UnitName ' Wire4']);

linkaxes([wire1, wire2, wire3, wire4], 'xy');

maxTemplate = eval(sprintf('wire%i', curUnitSummary.Spike_Features(1,1)));

text(maxTemplate, curUnitSummary.Spike_Features(2,1), curUnitSummary.Spike_Features(2,2),...
    {'\uparrow', '\bf \fontsize{12} Valley'}, 'horizontalalignment', 'center', 'verticalalignment', 'top');
text(maxTemplate, curUnitSummary.Spike_Features(3,1), curUnitSummary.Spike_Features(3,2),...
    {'\bf \fontsize{12} Peak', '\downarrow'}, 'horizontalalignment', 'center', 'verticalalignment', 'bottom');

%% Plot Autocorrelation
autoCorrBySpk = nan(length(curUniSpikeTimes), 160);
for spk = 1:length(curUniSpikeTimes)
    curRelativeSpikeTimes = curUniSpikeTimes(spk)-curUniSpikeTimes;
    curRelativeSpikeTimes(curRelativeSpikeTimes==0)=[];
    autoCorrBySpk(spk,:) = histcounts(curRelativeSpikeTimes, linspace(-1,1,161));
end
autoCorr = sum(autoCorrBySpk,1);

subplot(2,4,5:6)
bar(linspace(-1,1,160), autoCorr, 'black');
title({'Autocorrellogram', sprintf('Mean Spike Rate = %.02g spk/s', curUnitSummary.Mean_SpikeRate)});
axis tight
xlabel Time
ylabel Counts

%% Plot Text blah
subplot(2,4,7:8)
set(gca, 'ylim', [0 7], 'xlim', [0 20]);
axis off
text(0,5, '\bf\fontsize{10}Mean', 'horizontalalignment', 'left');
text(0,4, '\bf\fontsize{10}Median', 'horizontalalignment', 'left');
text(0,3, {'\bf\fontsize{10}StDev', '(Ang/Circ)'}, 'horizontalalignment', 'left');
text(0,2, {'\bf\fontsize{10}Vector', 'Length'}, 'horizontalalignment', 'left');
text(0,1, {'\bf\fontsize{10}R-Test', '(p/r-value)'}, 'horizontalalignment', 'left');

text(5,6, '\bf \fontsize{9} \theta', 'horizontalalignment', 'center');
text(5,5, ['\fontsize{8}' sprintf('%.02g', curUnitSummary.Spike_Phase_Relations.Theta.Mean)], 'horizontalalignment', 'center');
text(5,4, ['\fontsize{8}' sprintf('%.02g', curUnitSummary.Spike_Phase_Relations.Theta.Median)], 'horizontalalignment', 'center');
text(5,3, ['\fontsize{8}' sprintf('%.02g (%.02g)', curUnitSummary.Spike_Phase_Relations.Theta.AngDev,...
    curUnitSummary.Spike_Phase_Relations.Theta.CircStDev)], 'horizontalalignment', 'center');
text(5,2, ['\fontsize{8}' sprintf('%.02g', curUnitSummary.Spike_Phase_Relations.Theta.R_Length)], 'horizontalalignment', 'center');
text(5,1.2, ['\fontsize{8}' sprintf('r=%.02f', curUnitSummary.Spike_Phase_Relations.Theta.R_Test(2))], 'horizontalalignment', 'center');
text(5,0.8, ['\fontsize{8}' sprintf('(p=%.02f)', curUnitSummary.Spike_Phase_Relations.Theta.R_Test(1))], 'horizontalalignment', 'center');

text(8,6, '\bf \fontsize{9} Low \beta', 'horizontalalignment', 'center');
text(8,5, ['\fontsize{8}' sprintf('%.02g', curUnitSummary.Spike_Phase_Relations.LowBeta.Mean)], 'horizontalalignment', 'center');
text(8,4, ['\fontsize{8}' sprintf('%.02g', curUnitSummary.Spike_Phase_Relations.LowBeta.Median)], 'horizontalalignment', 'center');
text(8,3, ['\fontsize{8}' sprintf('%.02g (%.02g)', curUnitSummary.Spike_Phase_Relations.LowBeta.AngDev,...
    curUnitSummary.Spike_Phase_Relations.LowBeta.CircStDev)], 'horizontalalignment', 'center');
text(8,2, ['\fontsize{8}' sprintf('%.02g', curUnitSummary.Spike_Phase_Relations.LowBeta.R_Length)], 'horizontalalignment', 'center');
text(8,1.2, ['\fontsize{8}' sprintf('r=%.02f', curUnitSummary.Spike_Phase_Relations.LowBeta.R_Test(2))], 'horizontalalignment', 'center');
text(8,0.8, ['\fontsize{8}' sprintf('(p=%.02f)', curUnitSummary.Spike_Phase_Relations.LowBeta.R_Test(1))], 'horizontalalignment', 'center');

text(11,6, '\bf \fontsize{9} \beta', 'horizontalalignment', 'center');
text(11,5, ['\fontsize{8}' sprintf('%.02g', curUnitSummary.Spike_Phase_Relations.Beta.Mean)], 'horizontalalignment', 'center');
text(11,4, ['\fontsize{8}' sprintf('%.02g', curUnitSummary.Spike_Phase_Relations.Beta.Median)], 'horizontalalignment', 'center');
text(11,3, ['\fontsize{8}' sprintf('%.02g (%.02g)', curUnitSummary.Spike_Phase_Relations.Beta.AngDev,...
    curUnitSummary.Spike_Phase_Relations.Beta.CircStDev)], 'horizontalalignment', 'center');
text(11,2, ['\fontsize{8}' sprintf('%.02g', curUnitSummary.Spike_Phase_Relations.Beta.R_Length)], 'horizontalalignment', 'center');
text(11,1.2, ['\fontsize{8}' sprintf('r=%.02f', curUnitSummary.Spike_Phase_Relations.Beta.R_Test(2))], 'horizontalalignment', 'center');
text(11,0.8, ['\fontsize{8}' sprintf('(p=%.02f)', curUnitSummary.Spike_Phase_Relations.Beta.R_Test(1))], 'horizontalalignment', 'center');

text(14,6, '\bf \fontsize{9} Low \gamma', 'horizontalalignment', 'center');
text(14,5, ['\fontsize{8}' sprintf('%.02g', curUnitSummary.Spike_Phase_Relations.LowGamma.Mean)], 'horizontalalignment', 'center');
text(14,4, ['\fontsize{8}' sprintf('%.02g', curUnitSummary.Spike_Phase_Relations.LowGamma.Median)], 'horizontalalignment', 'center');
text(14,3, ['\fontsize{8}' sprintf('%.02g (%.02g)', curUnitSummary.Spike_Phase_Relations.LowGamma.AngDev,...
    curUnitSummary.Spike_Phase_Relations.LowGamma.CircStDev)], 'horizontalalignment', 'center');
text(14,2, ['\fontsize{8}' sprintf('%.02g', curUnitSummary.Spike_Phase_Relations.LowGamma.R_Length)], 'horizontalalignment', 'center');
text(14,1.2, ['\fontsize{8}' sprintf('r=%.02f', curUnitSummary.Spike_Phase_Relations.LowGamma.R_Test(2))], 'horizontalalignment', 'center');
text(14,0.8, ['\fontsize{8}' sprintf('(p=%.02f)', curUnitSummary.Spike_Phase_Relations.LowGamma.R_Test(1))], 'horizontalalignment', 'center');

text(17,6, '\bf \fontsize{9} High \gamma', 'horizontalalignment', 'center');
text(17,5, ['\fontsize{8}' sprintf('%.02g', curUnitSummary.Spike_Phase_Relations.HighGamma.Mean)], 'horizontalalignment', 'center');
text(17,4, ['\fontsize{8}' sprintf('%.02g', curUnitSummary.Spike_Phase_Relations.HighGamma.Median)], 'horizontalalignment', 'center');
text(17,3, ['\fontsize{8}' sprintf('%.02g (%.02g)', curUnitSummary.Spike_Phase_Relations.HighGamma.AngDev,...
    curUnitSummary.Spike_Phase_Relations.HighGamma.CircStDev)], 'horizontalalignment', 'center');
text(17,2, ['\fontsize{8}' sprintf('%.02g', curUnitSummary.Spike_Phase_Relations.HighGamma.R_Length)], 'horizontalalignment', 'center');
text(17,1.2, ['\fontsize{8}' sprintf('r=%.02f', curUnitSummary.Spike_Phase_Relations.HighGamma.R_Test(2))], 'horizontalalignment', 'center');
text(17,0.8, ['\fontsize{8}' sprintf('(p=%.02f)', curUnitSummary.Spike_Phase_Relations.HighGamma.R_Test(1))], 'horizontalalignment', 'center');

text(20,6, '\bf \fontsize{9} Ripple', 'horizontalalignment', 'center');
text(20,5, ['\fontsize{8}' sprintf('%.02g', curUnitSummary.Spike_Phase_Relations.Ripple.Mean)], 'horizontalalignment', 'center');
text(20,4, ['\fontsize{8}' sprintf('%.02g', curUnitSummary.Spike_Phase_Relations.Ripple.Median)], 'horizontalalignment', 'center');
text(20,3, ['\fontsize{8}' sprintf('%.02g (%.02g)', curUnitSummary.Spike_Phase_Relations.Ripple.AngDev,...
    curUnitSummary.Spike_Phase_Relations.Ripple.CircStDev)], 'horizontalalignment', 'center');
text(20,2, ['\fontsize{8}' sprintf('%.02g', curUnitSummary.Spike_Phase_Relations.Ripple.R_Length)], 'horizontalalignment', 'center');
text(20,1.2, ['\fontsize{8}' sprintf('r=%.02f', curUnitSummary.Spike_Phase_Relations.Ripple.R_Test(2))], 'horizontalalignment', 'center');
text(20,0.8, ['\fontsize{8}' sprintf('(p=%.02f)', curUnitSummary.Spike_Phase_Relations.Ripple.R_Test(1))], 'horizontalalignment', 'center');

title({'Session Wide','Spike Phase Relationship'})

if saveYN==1
    set(figID, 'PaperOrientation', 'landscape');
    print('-fillpage', figID, '-dpdf', [curUnitSummary.UnitName '_Unit_Feature_Summary.pdf']);
end
