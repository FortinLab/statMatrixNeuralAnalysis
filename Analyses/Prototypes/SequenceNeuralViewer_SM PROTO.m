function SequenceNeuralViewer_SM
%%

%% Load data
[smFile,smPath] = uigetfile('*.mat');
cd(smPath);
files = dir(smPath);
fileNames = {files.name};
% Load the behavior matrix file for poke events plots
behMatFile = fileNames{cellfun(@(a)~isempty(a), strfind(fileNames, 'BehaviorMatrix'))};
nsmblMatFile = fileNames{cellfun(@(a)~isempty(a), strfind(fileNames, 'EnsembleMatrix'))};
load([smPath behMatFile]);
load([smPath nsmblMatFile]);
% Identify list of all statMatrix files
smFileList = fileNames(cellfun(@(a)~isempty(a), regexp(fileNames, '_SM\>')))';

%% Organize Ensemble Data
% Declare histogram bin size for peak FR latency decision
histBins = 0:0.1:1.2;

% Extract session variables & ensemble FR
behavMatrixTrialStruct = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [min(histBins) max(histBins)], 'PokeIn');
ensembleTrialData = ExtractTrialData_SM(behavMatrixTrialStruct, ensembleMatrix(:,2:end)); %#ok<*NODEF>
trialTSvect = behavMatrix(behavMatrixTrialStruct(1).TrialLogVect,1)-behavMatrix(behavMatrixTrialStruct(1).PokeInIndex,1);
inSeqLog = [behavMatrixTrialStruct.TranspositionDistance]==0;
perfLog = [behavMatrixTrialStruct.Performance]==1;
oAlog = [behavMatrixTrialStruct.Position]==1;
% iscEnsembleData = ensembleTrialData(inSeqLog & perfLog);
iscEnsembleData = ensembleTrialData(inSeqLog & perfLog & ~oAlog);

% Bin trial wise FR
iscEnsembleDataBinned = nan(length(histBins)-1, length(ensembleUnitSummaries), length(iscEnsembleData));
for tr = 1:length(iscEnsembleData)
    for uni = 1:length(ensembleUnitSummaries)
        iscEnsembleDataBinned(:,uni,tr) = histcounts(trialTSvect(logical(iscEnsembleData{tr}(:,uni))), histBins, 'Normalization', 'CountDensity');
    end
end
meanISCbinned = mean(iscEnsembleDataBinned,3);

% Normalize the FR and identify the peak FR bin
iscBinnedNormed = nan(size(meanISCbinned));
iscPeakLat = nan(size(meanISCbinned,2),3);
for uni = 1:size(meanISCbinned,2);
    iscBinnedNormed(:,uni) = meanISCbinned(:,uni)./max(meanISCbinned(:,uni));
    curPeakLat = find(iscBinnedNormed(:,uni)==1,1,'first');
    if isempty(curPeakLat)
        iscPeakLat(uni,:) = [size(iscBinnedNormed,1), uni, max(meanISCbinned(:,uni))];
    else
        iscPeakLat(uni,:) = [curPeakLat, uni, max(meanISCbinned(:,uni))];
    end
end
% Select cells based on peak FR value
trialFRthreshLog = iscPeakLat(:,3)<1;
iscBinnedNormed(:,trialFRthreshLog) = [];
iscPeakLat(trialFRthreshLog,:) = [];

% Separate out principal cells vs interneurons using spike width...
%%%%% NOTE: This value is chosen by eye based on the data currently... a
%%%%% more principled method for choosing the value should be used in the
%%%%% future
spkWdth = [ensembleUnitSummaries(~trialFRthreshLog).Spike_Width]*40;
spkRt = [ensembleUnitSummaries(~trialFRthreshLog).Mean_SpikeRate];
% Width Based
% widthThresh = 0.45;
% pcLog = (spkWdth>widthThresh;
% Spike Rate Based
spkRtThresh = 5;
pcLog = spkRt<spkRtThresh;


iscBinnedNormedPCs = iscBinnedNormed(:,pcLog);
iscPeakLatPCs = iscPeakLat(pcLog,:);
iscPeakLatPCsSorted = sortrows(iscPeakLatPCs);
iscBinnedNormedICs = iscBinnedNormed(:,~pcLog);
iscPeakLatICs = iscPeakLat(~pcLog,:);
iscPeakLatICsSorted = sortrows(iscPeakLatICs);

% Create Latency Alignments latency alignments
pcScatData = cell(2,size(iscPeakLatPCsSorted,1));
for uni = 1:size(iscPeakLatPCs,1)
    spikeLogVect = logical(ensembleMatrix(:,iscPeakLatPCsSorted(uni,2)+1));
    pcScatData{1,uni} = ensembleMatrix(spikeLogVect,1);
    uniID = ones(size(ensembleMatrix,1),1)*uni;
%     uniID = ones(size(ensembleMatrix,1),1)*find(peakLatSorted(:,2)==uni);
    pcScatData{2,uni} = uniID(spikeLogVect);
end

icScatData = cell(2,size(iscPeakLatICs,1));
for uni = 1:size(iscPeakLatICsSorted,1)
    spikeLogVect = logical(ensembleMatrix(:,iscPeakLatICsSorted(uni,2)+1));
    icScatData{1,uni} = ensembleMatrix(spikeLogVect,1);
    uniID = ones(size(ensembleMatrix,1),1)*uni + size(iscPeakLatPCs,1)+5;
%     uniID = ones(size(ensembleMatrix,1),1)*find(peakLatSort(:,2)==uni) + size(iscPeakLatPCs,1)+5;
    icScatData{2,uni} = uniID(spikeLogVect);
end

uniScatData = [pcScatData, icScatData];


%%%%%% This section is useful for a simple visualization of peak aligned FR values.
% sorted=sortrows([iscPeakLatICs(:,1), iscBinnedNormedICs']);
% figure; 
% subplot(2,1,1);
% imagesc(histBins(2:end)-0.05, 1:size(sorted,1)-1, sorted(:,2:end));
% % set(gca, 'ydir', 'reverse');
% hold on; line([0 0], [0.5 size(sorted,1)-0.5], 'linewidth', 2, 'color', 'k');
% subplot(2,1,2);
% histogram(iscPeakLatICs(:,3),0:1:40);
%%%%%%
%% Plot the session data
band2plot = 'Beta';
% Create Figure
figure;

% Plot Reference LFP Values
load(smFileList{1})
lfpColIDs = find(cellfun(@(a)~isempty(a), regexp(statMatrixColIDs, '_LFP_[(A-Z)|(a-z)]*\>')));
lfpBands = cellfun(@(a,b)a(b:end), statMatrixColIDs(lfpColIDs), regexp(statMatrixColIDs(lfpColIDs), '[(A-Z)|(a-z)]*\>'), 'uniformoutput', 0);
lfpVals = statMatrix(:,lfpColIDs(strcmp(lfpBands, band2plot)));
lfpVals = lfpVals/(max(abs(lfpVals)));

plot(statMatrix(:,1), lfpVals-1, 'color', [0.4 0.4 0.4]);

% Plot Spikes
hold on;
for uni = 1:size(uniScatData,2)
    scatter(uniScatData{1,uni}, uniScatData{2,uni}, '*k');
end
set(gca, 'ydir', 'reverse', 'color', 'none');

% Patch Trial Periods
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
        'YData', [zeros(1,2), repmat(icScatData{2,end}(1)+1,[1,2])],...
        'FaceColor', patchColor, 'FaceAlpha', 0.15,...
        'EdgeColor', patchColor);
end

%%
figure;
scatter([ensembleUnitSummaries.Mean_SpikeRate], [ensembleUnitSummaries.Spike_Width]*40); set(gca, 'xscale', 'log')
%% Create Figure
seqViewFig = figure;
figAxes = axes(seqViewFig, 'position', [0.25, 0.1, 0.7, 0.85]);
set(figAxes, 'ydir', 'reverse', 'color', 'none');

% UI Control Buttons
selectFilebtn = uicontrol(artifactRemovalFig, 'Units', 'Normalized', 'Style', 'listbox', 'String', smFileList,...
    'Position', [0.0375,0.65,0.125,0.28],'Callback', @selectFile);

%%