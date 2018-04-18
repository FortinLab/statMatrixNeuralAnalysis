
close all
clear all
%%
origDir = cd;
[fileDir] = uigetdir(origDir);
if fileDir==0
    disp('Analysis Cancelled')
    return
else
    cd(fileDir)
end
dirContents = dir(fileDir);
fileNames = {dirContents.name};

%% Load Relevant Data
load(fileNames{cellfun(@(a)~isempty(a), strfind(fileNames, 'EnsembleMatrix'))});

%%
spr = [ensembleUnitSummaries.Spike_Phase_Relations];
figure;
theta = [spr.Theta]; theta = [theta.R_Length]; a=subplot(6,1,1); histogram(theta, 0:.01:.6); axis tight;
lowBeta = [spr.LowBeta]; lowBeta = [lowBeta.R_Length]; b=subplot(6,1,2); histogram(lowBeta, 0:.01:.6);axis tight;
beta = [spr.Beta]; beta = [beta.R_Length]; c=subplot(6,1,3); histogram(beta, 0:.01:.6);axis tight;
lowGamma = [spr.LowGamma]; lowGamma = [lowGamma.R_Length]; d=subplot(6,1,4); histogram(lowGamma, 0:.01:.6);axis tight;
highGamma = [spr.HighGamma]; highGamma = [highGamma.R_Length]; e=subplot(6,1,5); histogram(highGamma, 0:.01:.6);axis tight;
ripple = [spr.Ripple]; ripple = [ripple.R_Length]; f=subplot(6,1,6); histogram(ripple, 0:.01:.6);axis tight;

linkaxes([a b c d e f], 'xy')

figure;
subplot(6,6,1); corrScatPlot(theta, theta, 'Theta','Theta',[]);
subplot(6,6,2); corrScatPlot(theta, lowBeta, 'Theta','LowBeta',[]);
subplot(6,6,3); corrScatPlot(theta, beta, 'Theta','Beta',[]);
subplot(6,6,4); corrScatPlot(theta, lowGamma, 'Theta','LowGamma',[]);
subplot(6,6,5); corrScatPlot(theta, highGamma, 'Theta','HighGamma',[]);
subplot(6,6,6); corrScatPlot(theta, ripple, 'Theta','Ripple',[]);

subplot(6,6,7); corrScatPlot(lowBeta, theta, 'LowBeta','Theta',[]);
subplot(6,6,8); corrScatPlot(lowBeta, lowBeta, 'LowBeta','LowBeta',[]);
subplot(6,6,9); corrScatPlot(lowBeta, beta, 'LowBeta','Beta',[]);
subplot(6,6,10); corrScatPlot(lowBeta, lowGamma, 'LowBeta','LowGamma',[]);
subplot(6,6,11); corrScatPlot(lowBeta, highGamma, 'LowBeta','HighGamma',[]);
subplot(6,6,12); corrScatPlot(lowBeta, ripple, 'LowBeta','Ripple',[]);

subplot(6,6,13); corrScatPlot(beta, theta, 'Beta','Theta',[]);
subplot(6,6,14); corrScatPlot(beta, lowBeta, 'Beta','LowBeta',[]);
subplot(6,6,15); corrScatPlot(beta, beta, 'Beta','Beta',[]);
subplot(6,6,16); corrScatPlot(beta, lowGamma, 'Beta','LowGamma',[]);
subplot(6,6,17); corrScatPlot(beta, highGamma, 'Beta','HighGamma',[]);
subplot(6,6,18); corrScatPlot(beta, ripple, 'Beta','Ripple',[]);

subplot(6,6,19); corrScatPlot(lowGamma, theta, 'LowGamma','Theta',[]);
subplot(6,6,20); corrScatPlot(lowGamma, lowBeta, 'LowGamma','LowBeta',[]);
subplot(6,6,21); corrScatPlot(lowGamma, beta, 'LowGamma','Beta',[]);
subplot(6,6,22); corrScatPlot(lowGamma, lowGamma, 'LowGamma','LowGamma',[]);
subplot(6,6,23); corrScatPlot(lowGamma, highGamma, 'LowGamma','HighGamma',[]);
subplot(6,6,24); corrScatPlot(lowGamma, ripple, 'LowGamma','Ripple',[]);

subplot(6,6,25); corrScatPlot(highGamma, theta, 'HighGamma','Theta',[]);
subplot(6,6,26); corrScatPlot(highGamma, lowBeta, 'HighGamma','LowBeta',[]);
subplot(6,6,27); corrScatPlot(highGamma, beta, 'HighGamma','Beta',[]);
subplot(6,6,28); corrScatPlot(highGamma, lowGamma, 'HighGamma','LowGamma',[]);
subplot(6,6,29); corrScatPlot(highGamma, highGamma, 'HighGamma','HighGamma',[]);
subplot(6,6,30); corrScatPlot(highGamma, ripple, 'HighGamma','Ripple',[]);

subplot(6,6,31); corrScatPlot(ripple, theta, 'Ripple','Theta',[]);
subplot(6,6,32); corrScatPlot(ripple, lowBeta, 'Ripple','LowBeta',[]);
subplot(6,6,33); corrScatPlot(ripple, beta, 'Ripple','Beta',[]);
subplot(6,6,34); corrScatPlot(ripple, lowGamma, 'Ripple','LowGamma',[]);
subplot(6,6,35); corrScatPlot(ripple, highGamma, 'Ripple','HighGamma',[]);
subplot(6,6,36); corrScatPlot(ripple, ripple, 'Ripple','Ripple',[]);
