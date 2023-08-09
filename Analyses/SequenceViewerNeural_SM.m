% SequenceViewerNeural_SM
close all
clear all %#ok<CLALL>
clc
%% Initialize The Data For the Figure
global svnFig figAxes lfpPlotHandle lfpBand2Plot rasterHandles pCellList intCellList numSeqs seqWindows curSeqNumInput behavMatrixTrialStruct smFile
% Get the file
[smFile,smPath] = uigetfile('*.mat');
if smPath==0
    return;
end
% Change directory
cd(smPath);
% Load the file
load(smFile);
% Identify the LFP columns
lfpNdxNd = cellfun(@(a)regexp(a,'_LFP_', 'end'), statMatrixColIDs, 'uniformoutput', 0);
% Identify the hilbert transform columns...
hilbCols = cell2mat(cellfun(@(a)~isempty(a), strfind(statMatrixColIDs, 'HilbVals'), 'uniformoutput', 0));
lfpColLog = cellfun(@(a)~isempty(a), lfpNdxNd);
% ... so they can be removed from the selection
lfpColLog(hilbCols) = false;
lfpData = statMatrix(:,lfpColLog);
% Now identify the band names so they can be put in the list
lfpNdxNd(~lfpColLog) = [];
lfpColIDs = statMatrixColIDs(lfpColLog);
lfpBandNames = cellfun(@(a,b)a(b+1:end),lfpColIDs, lfpNdxNd, 'uniformoutput', 0)';

% Load the behavior and ensemble matrices
files = dir(smPath);
fileNames = {files.name};
behMatFile = fileNames{cellfun(@(a)~isempty(a), strfind(fileNames, 'BehaviorMatrix'))};
nsmblMatFile = fileNames{cellfun(@(a)~isempty(a), strfind(fileNames, 'EnsembleMatrix'))};
load([smPath behMatFile]);
load([smPath nsmblMatFile]);
% Extract spike times and IDs for each unit
unitSpikeTimes = cell(size(ensembleMatrix,2)-1,1);
for uni = 2:size(ensembleMatrix,2)
    unitSpikeTimes{uni-1} = ensembleMatrix(ensembleMatrix(:,uni)~=0,1);
end
unitIDs = ensembleMatrixColIDs(:,2:end)';

% To find the time window to use we're going to identify how long the rat
% stayed in the port on average for InSeq correct trials and use that
% window as our timewindow for alinging things to.
% First pull out the poke in indices
trlStruct = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [0 0], 'PokeIn');
% Pull out relevant logicals
seqLog = [trlStruct.TranspositionDistance] == 0;
perfLog = [trlStruct.Performance] == 1;
% Extract sequence beginning and end points with 1s padding on either side
numSeqs = length(unique([trlStruct.SequenceNum]));
seqWindows = nan(numSeqs,2);
for seq = 1:numSeqs
    curSeq = trlStruct([trlStruct.SequenceNum]==seq);
    seqWindows(seq,1) = behavMatrix(curSeq(1).PokeInIndex-1500,1);
    seqWindows(seq,2) = behavMatrix(curSeq(end).PokeOutIndex+1500,1);
end
% Now determine the mean poke duration
pokeDur = [trlStruct.PokeOutIndex] - [trlStruct.PokeInIndex];
meanPokeDurISC = mean(pokeDur(seqLog & perfLog))/1000;
trialWindows = [-1.2 ceil(meanPokeDurISC*10)/10];
% Use the Organize/Extract functions to pull out and organize the trial data 
behavMatrixTrialStruct = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, trialWindows, 'PokeIn');
ensembleTrialData = ExtractTrialData_SM(behavMatrixTrialStruct, ensembleMatrix(:,2:end)); %#ok<*NODEF>

% Make a timestamp vector for use with the trial relative spiking data
trialTSvect = behavMatrix(behavMatrixTrialStruct(1).TrialLogVect,1)-behavMatrix(behavMatrixTrialStruct(1).PokeInIndex,1);
% Create a timestamp vector
histBins = trialWindows(1):0.1:trialWindows(2);
trialEnsembleDataBinned = nan(length(histBins)-1, length(unitIDs), length(ensembleTrialData));
for tr = 1:length(ensembleTrialData)
    for uni = 1:length(unitIDs)
        trialEnsembleDataBinned(:,uni,tr) = histcounts(trialTSvect(logical(ensembleTrialData{tr}(:,uni))), histBins, 'Normalization', 'CountDensity');
    end
end

plotData = struct('TrialEnsembleBinned', {trialEnsembleDataBinned},...
    'UnitIDs', {unitIDs},...
    'UnitSpikeTimes', {unitSpikeTimes},...
    'SessionTime', {statMatrix(:,1)},...
    'InSeqLog', {[trlStruct.TranspositionDistance]==0},...
    'TrialOdor', {[trlStruct.Odor]},...
    'TrialPosition', {[trlStruct.Position]},...
    'Performance', {[trlStruct.Performance]},...
    'UnitInfo', {ensembleUnitSummaries});
%% 
svnFig = figure;
figAxes = axes(svnFig, 'position', [0.3, 0.2, 0.65, 0.75]);
set(svnFig, 'userData', plotData);

ensembleData = DetermineEnsembleOrganization;

[pcData, intData] = SeparatePCandINTs;

%% Make the initial Figure
lfpVals = statMatrix(:,2);
lfpVals = (lfpVals/(max(abs(lfpVals))))*2;
lfpPlotHandle = plot(figAxes, statMatrix(:,1), lfpVals, 'color', [0.4 0.4 0.4]);
hold on;
rasterHandles = nan(length(unitIDs)+3,1);
for r = 1:length(rasterHandles)
    rasterHandles(r) = scatter(figAxes,1,(r+1)*-1, 'markerfacecolor', 'k', 'markeredgecolor', 'none');
end
UpdateRasters(pcData);
UpdateRasters(intData);

for trl = 1:length(behavMatrixTrialStruct)
    switch behavMatrixTrialStruct(trl).Odor
        case 1
            patchColor = [44/255 168/255 224/255];
        case 2
            patchColor = [154/255 133/255 122/255];
        case 3
            patchColor = [9/255 161/255 74/255];
        case 4 
            patchColor = [128/255 66/255 151/255];
        case 5
            patchColor = 'k';
    end
    patch(gca, 'XData', [behavMatrix(behavMatrixTrialStruct(trl).PokeInIndex,1), behavMatrix(behavMatrixTrialStruct(trl).PokeOutIndex,1), behavMatrix(behavMatrixTrialStruct(trl).PokeOutIndex,1),behavMatrix(behavMatrixTrialStruct(trl).PokeInIndex,1)],...
        'YData', [ones(1,2), repmat(length(plotData.UnitIDs)+1,[1,2])*-1],...
        'FaceColor', patchColor, 'FaceAlpha', 0.15,...
        'EdgeColor', patchColor);
end
refreshdata(svnFig);
   
axis tight
title(figAxes, smFile, 'interpreter', 'none');
set(figAxes, 'yticklabel', []);
set(figAxes, 'xlim', seqWindows(1,:));

%% Create the Rest of the Figure
% Choose LFP File Button
chooseLFPfileBtn = uicontrol(svnFig, 'Units', 'Normalized', 'Style', 'pushbutton', 'String', 'Choose LFP File',...
    'Position', [0.05,0.9,0.2,0.05], 'Tag', 'chooseLFPfile_Btn', 'Callback', @ChooseLFPfile);
% LFP Band Listbox
lfpBandList = axes(svnFig, 'position', [0.05, 0.875, 0.08, 0.05]);
set(lfpBandList, 'xlim', [-0.5 0.5], 'ylim', [0 0.5]);
text(lfpBandList, 0,0, 'LFP Data Options', 'horizontalalignment', 'center', 'verticalalignment', 'bottom');
axis(lfpBandList, 'off');
lfpBand2Plot = uicontrol(svnFig, 'Units', 'Normalized', 'Style', 'listbox', 'String', lfpBandNames,...
    'Position', [0.05,0.6,0.08,0.275], 'Tag', 'bandToPlot', 'Callback', @ChooseLFPbandToPlot, 'userData', lfpData);

% Temporal Organization Button Group
tempOrgMethod = uibuttongroup(svnFig, 'Units', 'Normalized', 'Position', [0.05, 0.43, 0.085, 0.15],...
    'Title', 'Organization Method', 'TitlePosition', 'centertop', 'Tag', 'orgMthd_BtnGrp',...
    'SelectionChangedFcn', @ChangeTempOrgMethod);
% Temp Org by Mean A
meanAorg = uicontrol(tempOrgMethod, 'Units', 'Normalized', 'Style', 'radiobutton', 'String', 'Peak Seq FR (reset with new seq)',...
    'Position', [0.05, 0.7, 0.9, 0.2], 'Tag', 'orgMeanA_PB');
% Temp Org by Mean InSeq
meanISorg = uicontrol(tempOrgMethod, 'Units', 'Normalized', 'Style', 'radiobutton', 'String', 'Mean InSeq',...
    'Position', [0.05, 0.4, 0.9, 0.2], 'Tag', 'orgMeanIS_PB');
% Temp Org by Non A InSeq
meanISsfpOrg = uicontrol(tempOrgMethod, 'Units', 'Normalized', 'Style', 'radiobutton', 'String', 'IS-SFP',...
    'Position', [0.05, 0.1, 0.9, 0.2], 'Tag', 'orgMeanISsfp_PB');

% Criteria for PC/IN Classification
cellTypeWidthClassify = uicontrol(svnFig, 'Units', 'Normalized', 'Style', 'pushbutton', 'String', 'Spike Width Org',...
    'Position', [0.05, 0.35, 0.052, 0.06], 'Tag', 'classCellWidth', 'Callback', @SpikeWidthOrganize);
cellTypeWidthThresh = uicontrol(svnFig, 'Units', 'Normalized', 'Style', 'edit', 'String', '0.6',...
    'Position', [0.105, 0.35, 0.025, 0.06]);

cellTypeRateClassify = uicontrol(svnFig, 'Units', 'Normalized', 'Style', 'pushbutton', 'String', 'Spike Rate Org',...
    'Position', [0.05 0.29, 0.052, 0.06], 'Tag', 'classCellRate', 'Callback', @SpikeRateOrganize);
cellTypeRateThresh = uicontrol(svnFig, 'Units', 'Normalized', 'Style', 'edit', 'String', '5',...
    'Position', [0.105, 0.29, 0.025, 0.06]);

% Principal Cells List
pcListTitleAx = axes(svnFig, 'position', [0.14, 0.875, 0.11, 0.05]);
set(pcListTitleAx, 'xlim', [-0.5 0.5], 'ylim', [0 0.5]);
text(pcListTitleAx, 0,0, 'Principal Cells', 'horizontalalignment', 'center', 'verticalalignment', 'bottom');
axis(pcListTitleAx, 'off');
pCellList = uicontrol(svnFig, 'Units', 'Normalized', 'Style', 'listbox', 'String', pcData(:,3),...
    'Position', [0.14, 0.575, 0.11, 0.3], 'Tag', 'pcCells_Lst', 'Callback', @SelectPC, 'userData', pcData);

% Move PC to Int list
movePCtoInt = uicontrol(svnFig, 'Units', 'Normalized', 'Style', 'pushbutton', 'String', 'PC>>Int',...
    'Position', [0.14, 0.5175, 0.05, 0.05], 'Callback', @PCtoINT);
% Move Int to PC list
moveInttoPC = uicontrol(svnFig, 'Units', 'Normalized', 'Style', 'pushbutton', 'String', 'PC<<Int',...
    'Position', [0.2, 0.5175, 0.05, 0.05], 'Callback', @INTtoPC);

% Interneurons List
intListTitleAx = axes(svnFig, 'position', [0.14, 0.5, 0.11, 0.05]);
set(intListTitleAx, 'xlim', [-0.5 0.5], 'ylim', [0 0.5]);
text(intListTitleAx, 0,0, 'Interneurons', 'horizontalalignment', 'center', 'verticalalignment', 'bottom');
axis(intListTitleAx, 'off');
intCellList = uicontrol(svnFig, 'Units', 'Normalized', 'Style', 'listbox', 'String', intData(:,3),...
    'Position', [0.14, 0.2, 0.11, 0.3], 'Tag', 'intCells_Lst', 'Callback', @SelectInt, 'userData', intData);

% Remove PC Button
rmvPC = uicontrol(svnFig, 'Units', 'Normalized', 'Style', 'pushbutton', 'String', 'Remove PC',...
    'Position', [0.05, 0.25, 0.08, 0.025], 'Callback', @RemovePC);

% Remove Int Button
rmvINT = uicontrol(svnFig, 'Units', 'Normalized', 'Style', 'pushbutton', 'String', 'Remove INT',...
    'Position', [0.05, 0.2, 0.08, 0.025], 'Callback', @RemoveINT);

% Popout Button
popOutButton = uicontrol(svnFig, 'Units', 'Normalized', 'Style', 'pushbutton', 'String', 'Pop-Out Figure',...
    'Position', [0.05, 0.05, 0.2, 0.1], 'Callback', @PopoutPlot);

% Previous Sequence Button
prevSeqBtn = uicontrol(svnFig, 'Units', 'Normalized', 'Style', 'pushbutton', 'String', '<< Previous Sequence',...
    'Position', [0.3, 0.05, 0.2, 0.1], 'Callback', @PrevSeq);
% Next Sequence Button
nextSeqBtn = uicontrol(svnFig, 'Units', 'Normalized', 'Style', 'pushbutton', 'String', 'Next Sequence >>',...
    'Position', [0.75, 0.05, 0.2, 0.1], 'Callback', @NextSeq);

% Sequence Number
curSeqNumInput = uicontrol(svnFig, 'Units', 'Normalized', 'Style', 'edit', 'String', '1',...
    'Position', [0.575, 0.075, 0.05, 0.05], 'Callback', @SelectSeq);
totalSeqNumAxs = axes('position', [0.625, 0.075, 0.05, 0.05]);
set(totalSeqNumAxs, 'xlim', [0 0.5], 'ylim', [-0.5 0.5]);
seqText = text(totalSeqNumAxs, 0,0, sprintf('/%i', numSeqs), 'horizontalalignment', 'left', 'verticalalignment', 'middle');
axis(totalSeqNumAxs, 'off');

%%
function [spikeData] = DetermineEnsembleOrganization(orgType)
global svnFig figAxes
seqWin = round(get(figAxes, 'xlim'));
plotData = get(svnFig, 'userData');
spikeTimes = plotData.UnitSpikeTimes;
unitIDs = plotData.UnitIDs;
unitInfo = arrayfun(@(a){a}, plotData.UnitInfo);
if nargin==0
    orgType = 1;
end
switch orgType
    case 1
        trialLog = plotData.TrialOdor==1 & plotData.TrialPosition==1 & plotData.Performance==1;
    case 2
        trialLog = plotData.InSeqLog==1 & plotData.Performance==1;
    case 3
        trialLog = plotData.InSeqLog==1 & plotData.Performance==1 & plotData.TrialOdor~=1;
end
if orgType == 2 || orgType == 3
    trialHist = mean(plotData.TrialEnsembleBinned(:,:,trialLog),3)';
    normHist = trialHist./max(trialHist,[],2);
else
    instFRgauss = (gausswin(200))/0.2;
    normHist = cell2mat(cellfun(@(a){conv(histcounts(a(a>min(seqWin) & a<max(seqWin)), plotData.SessionTime(plotData.SessionTime>min(seqWin) & plotData.SessionTime<max(seqWin))), instFRgauss, 'same')},spikeTimes));
    for u = 1:size(normHist,1)
        normHist(u,:) = normHist(u,:)./max(normHist(u,:));
    end
end
silentCellLog = sum(isnan(normHist),2)==size(normHist,2);
% lowFRlog = sum(trialHist<=1,2)==size(trialHist,2); silentCellLog = silentCellLog | lowFRlog;
spikeTimes(silentCellLog) = [];
unitIDs(silentCellLog) = [];
normHist(silentCellLog,:) = [];
unitInfo(:,silentCellLog) = [];
latency = cell(sum(~silentCellLog),1);
for u = 1:size(normHist,1)
    latency{u} = find(normHist(u,:)==1,1,'first');
end
spikeData = [latency, unitIDs, spikeTimes, unitInfo'];
plotData.SpikeData = spikeData;
set(svnFig, 'userData', plotData);
end
%%
function [pcData, intData] = SeparatePCandINTs(orgType, orgVal)
global svnFig
plotData = get(svnFig, 'userData');
spikeData = plotData.SpikeData;
uniFR = cellfun(@(a)a.Mean_SpikeRate, spikeData(:,4));
uniSW = cellfun(@(a)a.Spike_Width, spikeData(:,4))*40;
if nargin==0
    orgType = 'Width';
    orgVal = 0.6;
end

switch orgType
    case 'Width'
        pcLog = uniSW>=orgVal;
    case 'Rate'
        pcLog = uniFR<=orgVal;
end

% Select and Organize Putative Principal Cells
pcData = spikeData(pcLog,:);
pcSortVect = sortrows([cell2mat(pcData(:,1)), [1:sum(pcLog)]']); %#ok<NBRAK>
pcLines = 2:size(pcData,1)+1;
pcData = [num2cell(pcLines'), pcData(pcSortVect(:,2),:)];

% Select and Organize Putative Interneurons
intData = spikeData(~pcLog,:);
intSortVect = sortrows([cell2mat(intData(:,1)), [1:sum(~pcLog)]']); %#ok<NBRAK>
intLines = pcLines(end)+3:pcLines(end)+2+size(intData,1);
intData = [num2cell(intLines'), intData(intSortVect(:,2),:)];
end
%%
function UpdateRasters(dataMtx)
global rasterHandles
for u = 1:size(dataMtx,1)
    set(rasterHandles(dataMtx{u,1}), 'Xdata', dataMtx{u,4}, 'Ydata', ones(1,length(dataMtx{u,4})).*(dataMtx{u,1}*-1));
end
end
%%
function ChooseLFPfile(source,event)
global svnFig figAxes lfpPlotHandle lfpBand2Plot rasterHandles pCellList intCellList numSeqs seqWindows curSeqNumInput
[smFile,smPath] = uigetfile('*.mat');
if smPath==0
    return;
end

if strcmp(smPath, [cd '\'])
    load(smFile);
    lfpNdxNd = cellfun(@(a)regexp(a,'_LFP_', 'end'), statMatrixColIDs, 'uniformoutput', 0);
    hilbCols = cell2mat(cellfun(@(a)~isempty(a), strfind(statMatrixColIDs, 'HilbVals'), 'uniformoutput', 0));
    lfpColLog = cellfun(@(a)~isempty(a), lfpNdxNd);
    lfpColLog(hilbCols) = false;
    lfpData = statMatrix(:,lfpColLog);
    lfpNdxNd(~lfpColLog) = [];
    lfpColIDs = statMatrixColIDs(lfpColLog);
    lfpBandNames = cellfun(@(a,b)a(b+1:end),lfpColIDs, lfpNdxNd, 'uniformoutput', 0)';
    set(lfpBand2Plot, 'Value', 1);
    set(lfpBand2Plot, 'String', lfpBandNames);
    set(lfpBand2Plot, 'userData', lfpData);
    ChooseLFPbandToPlot
    title(figAxes, smFile);
end

end

%%
function ChooseLFPbandToPlot(source,event)
global lfpPlotHandle lfpBand2Plot
lfpData = get(lfpBand2Plot, 'userData');
band2Plot = get(lfpBand2Plot, 'value');
lfpData2Plot = lfpData(:,band2Plot);
lfpData2Plot = (lfpData2Plot/(max(abs(lfpData2Plot))))*2;
% set(lfpPlotHandle, 'Ydata', lfpData2Plot*-1);
set(lfpPlotHandle, 'Ydata', lfpData2Plot);
end

%%
function SpikeWidthOrganize(source,event)
error('Make me!')
end
%% 
function SpikeRateOrganize(source,event)
error('Make me!')
end
%%
function SelectPC(source,event)
global rasterHandles pCellList
pData = get(pCellList, 'userData');
curListPos = get(pCellList, 'Value');
if curListPos>size(pData,1)
    curListPos = size(pData,1);
end
set(rasterHandles(cell2mat(pData(:,1))), 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'k');
set(rasterHandles(pData{curListPos,1}), 'MarkerFaceColor', 'r');
end
%%
function RemovePC(source,event)
global rasterHandles pCellList intCellList
pData = get(pCellList, 'userData');
pStart = pData{1,1};
iData = get(intCellList, 'userData');
curListPos = get(pCellList, 'Value');

% Remove Index
pData(curListPos,:) = [];
% Resort PCs
pcSortMtx = sortrows([cell2mat(pData(:,2)), [1:size(pData,1)]']); %#ok<NBRAK>
pData = pData(pcSortMtx(:,2),:);
pData(:,1) = num2cell(pStart:size(pData,1)+1);
% Resort Ints
intSortMtx = sortrows([cell2mat(iData(:,2)), [1:size(iData,1)]']); %#ok<NBRAK>
iData = iData(intSortMtx(:,2),:);
iData(:,1) = num2cell(pData{end,1}+3:pData{end,1}+2+size(iData,1));

UpdateRasters(pData);
UpdateRasters(iData);
set(pCellList, 'userData', pData);
set(intCellList, 'userData', iData);
set(pCellList, 'String', pData(:,3));
if curListPos<=size(pData,1)
    set(pCellList, 'Value', curListPos-1);
elseif curListPos > size(pData,1)
    set(pCellList, 'Value', size(pData,1));
end

set(rasterHandles(~ismember(1:length(rasterHandles), cell2mat([pData(:,1);iData(:,1)]))), 'Xdata', 1, 'Ydata', 1);
SelectPC
end
%%
function SelectInt(source,event)
global rasterHandles intCellList
intData = get(intCellList, 'userData');
curListPos = get(intCellList, 'Value');
if curListPos>size(intData,1) || isempty(curListPos)
    curListPos = size(intData,1);
end

set(rasterHandles(cell2mat(intData(:,1))), 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'k');
set(rasterHandles(intData{curListPos,1}), 'MarkerFaceColor', [0.1 0.8 0.1]);
end
%%
function RemoveINT(source,event)
global intCellList pCellList rasterHandles
pData = get(pCellList, 'userData');
pStart = pData{1,1};
iData = get(intCellList, 'userData');
curListPos = get(intCellList, 'Value');

iData(curListPos,:) = [];
% Resort PCs
pcSortMtx = sortrows([cell2mat(pData(:,2)), [1:size(pData,1)]']); %#ok<NBRAK>
pData = pData(pcSortMtx(:,2),:);
pData(:,1) = num2cell(pStart:size(pData,1)+1);
% Resort Ints
intSortMtx = sortrows([cell2mat(iData(:,2)), [1:size(iData,1)]']); %#ok<NBRAK>
iData = iData(intSortMtx(:,2),:);
iData(:,1) = num2cell(pData{end,1}+3:pData{end,1}+2+size(iData,1));
UpdateRasters(iData);
UpdateRasters(pData);
set(pCellList, 'userData', pData);
set(intCellList, 'userData', iData);
set(intCellList, 'String', iData(:,3));
if curListPos<=length(pData)
    set(intCellList, 'Value', curListPos-1);
elseif curListPos > size(pData,1)
    set(intCellList, 'Value', size(pData,1));
end
set(rasterHandles(~ismember(1:length(rasterHandles), cell2mat([pData(:,1);iData(:,1)]))), 'Xdata', 1, 'Ydata', 1);
SelectInt
end
%%
function PCtoINT(source, event)
global rasterHandles pCellList intCellList
pData = get(pCellList, 'userData');
pStart = pData{1,1};
curPC = get(pCellList, 'Value');
pStr = get(pCellList, 'String');
cell2MoveStr = pStr(curPC);
iData = get(intCellList, 'userData');


% Modify Lists
iData = [iData;pData(curPC,:)];
pData(curPC,:) = [];

% Reorganize PCs
pcSortMtx = sortrows([cell2mat(pData(:,2)), [1:size(pData,1)]']); %#ok<NBRAK>
pData = pData(pcSortMtx(:,2),:);
pData(:,1) = num2cell(pStart:size(pData,1)+1);
for p = 1:size(pData,1)
    set(rasterHandles(pData{p,1}), 'Xdata', pData{p,4}, 'Ydata', ones(1,length(pData{p,4})).*(pData{p,1}*-1));
end
set(pCellList, 'String', pData(:,3));

% Reorganize INTs
intSortMtx = sortrows([cell2mat(iData(:,2)), [1:size(iData,1)]']); %#ok<NBRAK>
iData = iData(intSortMtx(:,2),:);
iData(:,1) = num2cell(pData{end,1}+3:pData{end,1}+2+size(iData,1));
for i = 1:size(iData,1)
    set(rasterHandles(iData{i,1}), 'Xdata', iData{i,4}, 'Ydata', ones(1,length(iData{i,4})).*(iData{i,1}*-1));
end
set(intCellList, 'String', iData(:,3));

set(pCellList, 'userData', pData);
set(intCellList, 'userData', iData);
set(intCellList, 'Value', find(strcmp(iData(:,3), cell2MoveStr)));
if curPC > 1
    set(pCellList, 'Value', curPC-1);
end

nonDataLog = ~ismember(1:length(rasterHandles), cell2mat([pData(:,1); iData(:,1)]));
set(rasterHandles(nonDataLog), 'Xdata', 1, 'Ydata', 1);
SelectInt
SelectPC
end
%%
function INTtoPC(source, event)
global rasterHandles pCellList intCellList
pData = get(pCellList, 'userData');
iData = get(intCellList, 'userData');
curInt = get(intCellList, 'Value');
iStr = get(intCellList, 'String');
cell2MoveStr = iStr(curInt);

% Modify Lists
pData = [pData;iData(curInt,:)];
ndxStart = pData{1,1};
iData(curInt,:) = [];

% Reorganize PCs
pcSortMtx = sortrows([cell2mat(pData(:,2)), [1:size(pData,1)]']); %#ok<NBRAK>
pData = pData(pcSortMtx(:,2),:);
pData(:,1) = num2cell(ndxStart:size(pData,1)+1);
for p = 1:size(pData,1)
    set(rasterHandles(pData{p,1}), 'Xdata', pData{p,4}, 'Ydata', ones(1,length(pData{p,4})).*(pData{p,1}*-1));
end

% Reorganize INTs
intSortMtx = sortrows([cell2mat(iData(:,2)), [1:size(iData,1)]']); %#ok<NBRAK>
iData = iData(intSortMtx(:,2),:);
iData(:,1) = num2cell(pData{end,1}+3:pData{end,1}+2+size(iData,1));
for i = 1:size(iData,1)
    set(rasterHandles(iData{i,1}), 'Xdata', iData{i,4}, 'Ydata', ones(1,length(iData{i,4})).*(iData{i,1}*-1));
end

set(pCellList, 'userData', pData);
set(pCellList, 'String', pData(:,3));
set(intCellList, 'userData', iData);
set(intCellList, 'String', iData(:,3));
set(pCellList, 'Value', find(strcmp(pData(:,3), cell2MoveStr)));
if curInt > 1
    set(intCellList, 'Value', curInt-1);
end

nonDataLog = ~ismember(1:length(rasterHandles), cell2mat([pData(:,1); iData(:,1)]));
set(rasterHandles(nonDataLog), 'Xdata', 1, 'Ydata', 1);
SelectPC
SelectInt
end
%%
function ChangeTempOrgMethod(source,event)
global figAxes pCellList intCellList rasterHandles
orgMethodVal = logical(cell2mat(get(get(source, 'Children'), 'Value')));
orgMethodTag = get(get(source, 'Children'), 'Tag');
curMethod = orgMethodTag{orgMethodVal};
switch curMethod
    case 'orgMeanA_PB'
        [spikeData] = DetermineEnsembleOrganization(1);
    case 'orgMeanIS_PB'
        [spikeData] = DetermineEnsembleOrganization(2);
    case 'orgMeanISsfp_PB'
        [spikeData] = DetermineEnsembleOrganization(3);
end
iData = get(intCellList, 'userData');
intRmvLog = ismember(spikeData(:,2), iData(:,3));
pData = spikeData;
pData(intRmvLog,:) = [];
pSortVect = sortrows([cell2mat(pData(:,1)), [1:size(pData,1)]']); %#ok<NBRAK>
pLines = 2:size(pData,1)+1;
pData = [num2cell(pLines'), pData(pSortVect(:,2),:)];
UpdateRasters(pData);
set(pCellList, 'userData', pData);
set(pCellList, 'String', pData(:,3));

for i = 1:size(iData,1)
    newSpikeData = spikeData(strcmp(iData{i,3}, spikeData(:,2)),:);
    if isempty(newSpikeData)
        iData{i,1} = nan;
    else
        iData(i,2:end) = newSpikeData;
    end
end
nanCheck = isnan(cell2mat(iData(:,1)));
iData(nanCheck,:) = [];
iData = iData(:,2:end);
iSortVect = sortrows([cell2mat(iData(:,1)), [1:size(iData,1)]']); %#ok<NBRAK>
iLines = pLines(end)+3:pLines(end)+2+size(iData,1);
iData = [num2cell(iLines'), iData(iSortVect(:,2),:)];
UpdateRasters(iData);
set(figAxes, 'ylim', [(iData{end,1}+1)*-1 1]);
set(intCellList, 'userData', iData);
set(intCellList, 'String', iData(:,3));
SelectPC
SelectInt

set(rasterHandles(~ismember(1:length(rasterHandles), cell2mat([pData(:,1);iData(:,1)]))), 'Xdata', 1, 'Ydata', 1);
end

%%
function PrevSeq(source,event)
global seqWindows curSeqNumInput figAxes
curSeqNum = str2double(get(curSeqNumInput, 'String'));
if curSeqNum - 1 >= 1
    set(curSeqNumInput, 'String', num2str(curSeqNum-1));
    set(figAxes, 'xlim', seqWindows(curSeqNum-1,:));
end
end
%%
function NextSeq(source,event)
global numSeqs seqWindows curSeqNumInput figAxes
curSeqNum = str2double(get(curSeqNumInput, 'String'));
if curSeqNum + 1 <= numSeqs
    set(curSeqNumInput, 'String', num2str(curSeqNum+1));
    set(figAxes, 'xlim', seqWindows(curSeqNum+1,:));
end
end
%%
function SelectSeq(source,event)
global numSeqs seqWindows curSeqNumInput figAxes
curSeqNum = str2double(get(curSeqNumInput, 'String'));
if curSeqNum <= numSeqs
    set(curSeqNumInput, 'String', num2str(curSeqNum));
    set(figAxes, 'xlim', seqWindows(curSeqNum,:));
end
end
%%
function PopoutPlot(source,event)
global svnFig lfpBand2Plot pCellList intCellList seqWindows curSeqNumInput behavMatrixTrialStruct smFile
plotData = get(svnFig, 'userData');
poFig = figure;
lfpData = get(lfpBand2Plot, 'userData');
curLFP = lfpData(:,get(lfpBand2Plot, 'Value'));
curLFP = (curLFP/(max(abs(curLFP))))*2;
plot(plotData.SessionTime, curLFP, 'color', [0.4 0.4 0.4]);
hold on;
pData = get(pCellList, 'userData');
iData = get(intCellList, 'userData');

scatterData = [pData(:,1), pData(:,4); iData(:,1), iData(:,4)];
for u = 1:size(scatterData,1)
    scatter(scatterData{u,2}, ones(size(scatterData{u,2})).*scatterData{u,1}*-1, 'markerfacecolor', 'k', 'markeredgecolor', 'none');
end

for trl = 1:length(behavMatrixTrialStruct)
    switch behavMatrixTrialStruct(trl).Odor
        case 1
            patchColor = [44/255 168/255 224/255];
        case 2
            patchColor = [154/255 133/255 122/255];
        case 3
            patchColor = [9/255 161/255 74/255];
        case 4 
            patchColor = [128/255 66/255 151/255];
        case 5
            patchColor = 'k';
    end
    patch(gca, 'XData', [plotData.SessionTime(behavMatrixTrialStruct(trl).PokeInIndex,1), plotData.SessionTime(behavMatrixTrialStruct(trl).PokeOutIndex,1), plotData.SessionTime(behavMatrixTrialStruct(trl).PokeOutIndex,1),plotData.SessionTime(behavMatrixTrialStruct(trl).PokeInIndex,1)],...
        'YData', [ones(1,2), repmat(iData{end,1}+1,[1,2])*-1],...
        'FaceColor', patchColor, 'FaceAlpha', 0.15,...
        'EdgeColor', patchColor);
end
set(gca, 'xlim', seqWindows(str2double(get(curSeqNumInput, 'String')),:), 'ylim', [-1*iData{end,1}-1, 2]);
title(sprintf('%s; Sequence %s', smFile, get(curSeqNumInput, 'String')), 'interpreter', 'none');
orient(gcf, 'tall');
orient(gcf, 'landscape');
set(gcf, 'Renderer', 'Painters');
end