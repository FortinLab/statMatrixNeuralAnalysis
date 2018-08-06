function [ensembleOutput] = EnsembleSummary_SM(fileDir)
%% EnsembleSummary PROTO
%
%
%%
if nargin==0
    [fileDir] = uigetdir;
    if fileDir==0
        disp('Analysis Cancelled')
        return
    else
        cd(fileDir)
    end
else
    cd(fileDir);
end
%%
dirContents = dir(fileDir);
fileNames = {dirContents.name};

load(fileNames{cellfun(@(a)~isempty(a), strfind(fileNames, 'BehaviorMatrix'))});
load(fileNames{cellfun(@(a)~isempty(a), strfind(fileNames, 'EnsembleMatrix'))});
ensembleUniDta = ensembleMatrix(:,2:end);
samp = mode(diff(ensembleMatrix(:,1)));
%% Spike Phase Relations
uniNames = {ensembleUnitSummaries.UnitName};
spkWdth = [ensembleUnitSummaries.Spike_Width];
spkRt = [ensembleUnitSummaries.Mean_SpikeRate];

spikePR = [ensembleUnitSummaries.Spike_Phase_Relations];
thetaPR = [spikePR.Theta];
lowBetaPR = [spikePR.LowBeta];
betaPR = [spikePR.Beta];
lowGammaPR = [spikePR.LowGamma];
highGammaPR = [spikePR.HighGamma];
ripplePR = [spikePR.Ripple];

%%
figure;
corrScatPlot(spkWdth, spkRt, 'Spike Width', 'Overall Spike Rate (spk/s)', []);
%%
figure;
subplot(6,3,1)
corrScatPlot(spkWdth, [thetaPR.Mean], 'Spike Width', '\theta Phase', []);
subplot(6,3,2)
corrScatPlot(spkWdth, [thetaPR.R_Length], 'Spike Width', 'Vector Length', []);
subplot(6,3,3)
corrScatPlot(spkWdth, cell2mat(cellfun(@(a)a(2), {thetaPR.R_Test}, 'uniformoutput', 0)), 'Spike Width', 'R-test Val', []);

subplot(6,3,4)
corrScatPlot(spkWdth, [lowBetaPR.Mean], 'Spike Width', 'Low \beta Phase', []);
subplot(6,3,5)
corrScatPlot(spkWdth, [lowBetaPR.R_Length], 'Spike Width', 'Vector Length', []);
subplot(6,3,6)
corrScatPlot(spkWdth, cell2mat(cellfun(@(a)a(2), {lowBetaPR.R_Test}, 'uniformoutput', 0)), 'Spike Width', 'R-test Val', []);

subplot(6,3,7)
corrScatPlot(spkWdth, [betaPR.Mean], 'Spike Width', '\beta Phase', []);
subplot(6,3,8)
corrScatPlot(spkWdth, [betaPR.R_Length], 'Spike Width', 'Vector Length', []);
subplot(6,3,9)
corrScatPlot(spkWdth, cell2mat(cellfun(@(a)a(2), {betaPR.R_Test}, 'uniformoutput', 0)), 'Spike Width', 'R-test Val', []);

subplot(6,3,10)
corrScatPlot(spkWdth, [lowGammaPR.Mean], 'Spike Width', 'Low \gamma Phase', []);
subplot(6,3,11)
corrScatPlot(spkWdth, [lowGammaPR.R_Length], 'Spike Width', 'Vector Length', []);
subplot(6,3,12)
corrScatPlot(spkWdth, cell2mat(cellfun(@(a)a(2), {lowGammaPR.R_Test}, 'uniformoutput', 0)), 'Spike Width', 'R-test Val', []);

subplot(6,3,13)
corrScatPlot(spkWdth, [highGammaPR.Mean], 'Spike Width', 'High \gamma Phase', []);
subplot(6,3,14)
corrScatPlot(spkWdth, [highGammaPR.R_Length], 'Spike Width', 'Vector Length', []);
subplot(6,3,15)
corrScatPlot(spkWdth, cell2mat(cellfun(@(a)a(2), {highGammaPR.R_Test}, 'uniformoutput', 0)), 'Spike Width', 'R-test Val', []);

subplot(6,3,16)
corrScatPlot(spkWdth, [ripplePR.Mean], 'Spike Width', 'Ripple Phase', []);
subplot(6,3,17)
corrScatPlot(spkWdth, [ripplePR.R_Length], 'Spike Width', 'Vector Length', []);
subplot(6,3,18)
corrScatPlot(spkWdth, cell2mat(cellfun(@(a)a(2), {ripplePR.R_Test}, 'uniformoutput', 0)), 'Spike Width', 'R-test Val', []);