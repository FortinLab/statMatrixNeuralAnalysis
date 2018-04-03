function CreateStatMatrixFromPLX_MS
%% CreateStatMatrixFromPLX_MS
% Description goes here.
% For use with Boston files, will need to create an Irvine version

%%
origCD = cd;
[origPlxFile, filePath] = uigetfile('.plx','Identify Original .PLX File');
if origPlxFile == 0
    disp('No file selected');
    return
end
outputFileName = inputdlg('Name Output File', 'origPlxFile', 1, {origPlxFile});
cd(filePath);

spikeFileDir = uigetdir(cd, 'Where are the cut MountainSort .PLX files?');
spikeFileDirFiles = dir([spikeFileDir '\']);
spkFiles = {spikeFileDirFiles.name};
spkFiles = spkFiles(logical(cellfun(@(a)~isempty(a), strfind(spkFiles, '.plx'))));

%% Behavioral Data
% [behaviorData] = SummarizePLXabbr_BOS(origPlxFile);
[behaviorData] = SummarizePLXabbr(origPlxFile);
[samp, ~, tetTS, fn, ~] = plx_ad_v(origPlxFile, 2);
tsVect = ((0:(fn(1)-1))/samp)+tetTS(1);
for fragNum = 2:length(tetTS)
    tsVect = [tsVect, ((0:(fn(fragNum)-1))/samp)+tetTS(fragNum)]; %#ok<AGROW>
end
seqLength = behaviorData.Summary.SequenceLength;
behVals = nan(length(tsVect), seqLength*2 + 5);
behDataHeaders = cell(1,seqLength*2 + 7);
ssnData = behaviorData.Raw;
for seq = 1:seqLength
    itemPresTimes = [ssnData([ssnData.SequenceItem]==seq).ItemPresentationTime];
    posPresTimes = [ssnData([ssnData.OrdinalPosition]==seq).ItemPresentationTime];
    inSeqPresTimes = intersect(itemPresTimes, posPresTimes);
    outSeqPresTimes = setdiff(itemPresTimes,posPresTimes);
    behVals(:,seq) = histcounts(itemPresTimes, [0 tsVect])';
    behVals(:,seq+seqLength) = histcounts(posPresTimes, [0 tsVect])';
    behDataHeaders{seq} = ['Odor' num2str(seq)];
    behDataHeaders{seq+seqLength} = ['Position' num2str(seq)];
    if seq == 1
        behVals(:,seqLength+seqLength+1) = histcounts(inSeqPresTimes, [0 tsVect])' - histcounts(outSeqPresTimes,[0 tsVect])';
    else
        behVals(:,seqLength+seqLength+1) = behVals(:,seqLength+seqLength+1) + histcounts(inSeqPresTimes, [0 tsVect])' - histcounts(outSeqPresTimes,[0 tsVect])';
    end
end
behDataHeaders{seq+seqLength+1} = 'InSeqLog';
itmPresTimes = [ssnData.ItemPresentationTime];
trialPerformance = [ssnData.Performance];
corrTrials = itmPresTimes(logical(trialPerformance));
corTrlHistCounts = histcounts(corrTrials, [0 tsVect])';
inCorrTrials = itmPresTimes(~logical(trialPerformance));
inCorTrlHistCounts = histcounts(inCorrTrials, [0 tsVect])';
behVals(:,seqLength+seqLength+2) = corTrlHistCounts + (inCorTrlHistCounts*-1);
behDataHeaders{seq+seqLength+2} = 'PerformanceLog';

behVals(:,seqLength+seqLength+3) = histcounts([ssnData.OdorTrigPokeTime], [0 tsVect])';
behVals(:,seqLength+seqLength+3) = behVals(:,seqLength+seqLength+3) - histcounts([ssnData.OdorPokeWithdrawTime], [0 tsVect])';
behDataHeaders{seq+seqLength+3} = 'PokeEvents';

behVals(:,seqLength+seqLength+4) = histcounts([ssnData.FrontRewardTime], [0 tsVect])';
behDataHeaders{seq+seqLength+4} = 'FrontReward';

behVals(:,seqLength+seqLength+5) = histcounts([ssnData.BackRewardTime], [0 tsVect])';
behDataHeaders{seq+seqLength+5} = 'BackReward';


[numChans, chanNames] = plx_event_names(origPlxFile);
findingStrobed = 1;
while findingStrobed
    for chan = 1:numChans
        curChan = chanNames(chan,:);
        intVals = double(curChan);
        valLim = find(~(intVals==0), 1, 'last');
        strobedChanLog = strcmp(curChan(1:valLim), 'Strobed');
        if strobedChanLog
            [~, strobedTS, strobedSV] = plx_event_ts(origPlxFile, curChan);
            [~, ~, ~, aniPosition] = plx_vt_interpret(strobedTS, strobedSV);
            aniPosHistBins = find(histcounts(aniPosition(:,1), [0 tsVect]));
            behVals(aniPosHistBins,seqLength+seqLength+6) = aniPosition(:,2);
            behVals(aniPosHistBins,seqLength+seqLength+7) = aniPosition(:,3);
            findingStrobed = 0;
        end
    end
end
behDataHeaders{seq+seqLength+6} = 'XvalRatMazePosition';
behDataHeaders{seq+seqLength+7} = 'YvalRatMazePosition';
behavMatrix = [tsVect', behVals]; %#ok<NASGU>
behavMatrixColIDs = [{'TimeBin'}, behDataHeaders]; %#ok<NASGU>
save([filePath outputFileName{1} '_BehaviorMatrix.mat'], 'behavMatrix', 'behavMatrixColIDs');
disp('Behavior data saved.');

%% Neural Data
[~, ~, ~, contCountFl] = plx_info(origPlxFile, 1);
[~, plxADchanNames] = plx_adchan_names(origPlxFile);
contChanLog = ~(contCountFl==0);
plxADchanNamesWcont = plxADchanNames(contChanLog',:);
chanNames = cell(size(plxADchanNamesWcont,1),1);
tetNames = cell(size(plxADchanNamesWcont,1),1);
chan2keep = false(size(chanNames));
for chan = 1:length(chanNames)
    [s, e] = regexp(plxADchanNamesWcont(chan,:), '([A-Z,a-z]*)([0-9]*)-([0-9]*)');
    if ~isempty(s)
        chan2keep(chan) = true;
        [ts,te] = regexp(plxADchanNamesWcont(chan,:), '([A-Z,a-z]*)([0-9]*)-');
        tetNames{chan} = plxADchanNamesWcont(chan,ts:te-1);
    end
    chanNames{chan} = plxADchanNamesWcont(chan,s:e);
end
tetChans = chanNames(chan2keep);
tetNames = unique(tetNames(chan2keep));
for tet = 1:length(tetNames)
    curTet = tetNames{tet};
    tetNumSpot = regexp(curTet, '[0-9]');
    tetNum = str2double(curTet(tetNumSpot:end));
    
    % Extract LFP Data
    adChan = tetChans{logical(cellfun(@(a)~isempty(a), strfind(tetChans,[curTet '-'])))};
    [samp, ~, tetTS, tetFN, tetV] = plx_ad_v(origPlxFile, adChan);
    tsVect = ((0:(tetFN(1)-1))/samp)+tetTS(1);
    for fragNum = 2:length(tetTS)
        tsVect = [tsVect, ((0:(tetFN(fragNum)-1))/samp)+tetTS(fragNum)]; %#ok<AGROW>
    end
    % Filter LFP Data
    [tetRawHilb, tetRawLFP] = PhaseFreqDetectAbbr(tetV,tsVect);
    [tetThetaHilb, tetThetaLFP] = PhaseFreqDetectAbbr(tetV,tsVect,4,12);
    [tetLowBetaHilb, tetLowBetaLFP] = PhaseFreqDetectAbbr(tetV,tsVect,13,19);
    [tetBetaHilb, tetBetaLFP] = PhaseFreqDetectAbbr(tetV,tsVect,20,40);
    [tetLowGammaHilb, tetLowGammaLFP] = PhaseFreqDetectAbbr(tetV,tsVect,41,59);
    [tetHighGammaHilb, tetHighGammaLFP] = PhaseFreqDetectAbbr(tetV,tsVect,61,80);
    [tetRipHilb, tetRipLFP] = PhaseFreqDetectAbbr(tetV,tsVect,150,250);
    
    % Identify Spike File (if they exist)
    curSpkFileLog = cellfun(@(a,b)~isempty(a) | ~isempty(b), strfind(spkFiles, sprintf('tet%i-', tetNum)), strfind(spkFiles, sprintf('tet%02i-', tetNum)));
    if sum(curSpkFileLog)==1
        curSpkFile = [spikeFileDir '\' spkFiles{curSpkFileLog}];
        [spkFileTScount, ~, ~, ~] = plx_info(curSpkFile, 1);
        curNumUnis = sum(spkFileTScount(2:end,2)>0);
    elseif sum(curSpkFileLog)==0
        curNumUnis = 0;
    elseif sum(curSpkFileLog)>=2
        error('Multiple files from the same tetrode, check the file directory')
    end
    tetUnitTSs = cell(1,curNumUnis);
    tetUnitSRs = cell(1,curNumUnis);
    tetUnitIDs = cell(1,curNumUnis);
    for uni = 1:curNumUnis
        [~, tempTSs] = plx_ts(curSpkFile, 1, uni);
        tetUnitTSs{uni} = histcounts(tempTSs, [0 tsVect])';
        tetUnitIDs{uni} = [curTet '-U' num2str(uni)];
    end
    statMatrix = [tsVect',tetRawLFP, tetRawHilb,...
        tetThetaLFP, tetThetaHilb,...
        tetLowBetaLFP, tetLowBetaHilb,...
        tetBetaLFP, tetBetaHilb,...          
        tetLowGammaLFP, tetLowGammaHilb,...
        tetHighGammaLFP, tetHighGammaHilb,...
        tetRipLFP, tetRipHilb,...
        cell2mat(tetUnitTSs), cell2mat(tetUnitSRs), behVals]; %#ok<NASGU>
    statMatrixColIDs = [{'TimeBin'}, {[curTet '_LFP_Raw']}, {[curTet '_LFP_Raw_HilbVals']},...
        {[curTet '_LFP_Theta']}, {[curTet '_LFP_Theta_HilbVals']},...
        {[curTet '_LFP_LowBeta']}, {[curTet '_LFP_LowBeta_HilbVals']},...
        {[curTet '_LFP_Beta']}, {[curTet '_LFP_Beta_HilbVals']},...
        {[curTet '_LFP_LowGamma']}, {[curTet '_LFP_LowGamma_HilbVals']},...
        {[curTet '_LFP_HighGamma']}, {[curTet '_LFP_HighGamma_HilbVals']},...
        {[curTet '_LFP_Ripple']}, {[curTet '_LFP_Ripple_HilbVals']},...
        tetUnitIDs, behDataHeaders]; %#ok<NASGU>
    save([filePath outputFileName{1} '_' curTet '.mat'], 'statMatrix', 'statMatrixColIDs');
    disp([tetNames{tet} ' saved!']);
end

disp('Complete')
cd(origCD);
