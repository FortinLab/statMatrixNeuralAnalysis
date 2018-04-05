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
cd(filePath);

spikeFileDir = uigetdir(cd, 'Where are the cut MountainSort .PLX files?');
spikeFileDirFiles = dir([spikeFileDir '\']);
spkFiles = {spikeFileDirFiles.name};
spkFiles = spkFiles(logical(cellfun(@(a)~isempty(a), strfind(spkFiles, '.plx'))));

outputFileName = inputdlg('Name Output File', 'origPlxFile', 1, [origPlxFile(1:end-4) '_MS']);

%% Behavioral Data
[behaviorData] = SummarizePLXabbr_BOS(origPlxFile);
% [behaviorData] = SummarizePLXabbr(origPlxFile);
[samp, ~, tetTS, fn, ~] = plx_ad_v(origPlxFile, 2);
tsVect = ((0:(fn(1)-1))/samp)+tetTS(1);
for fragNum = 2:length(tetTS)
    tsVect = [tsVect, ((0:(fn(fragNum)-1))/samp)+tetTS(fragNum)]; %#ok<AGROW>
end
seqLength = behaviorData.Summary.SequenceLength;
if isfield(behaviorData.Raw, 'ErrorSignalTime')
    behDataHeaders = cell(1,seqLength*2 + 8);
    behVals = nan(length(tsVect), seqLength*2 + 8);
else
    behDataHeaders = cell(1,seqLength*2 + 7);
    behVals = nan(length(tsVect), seqLength*2 + 7);
end
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
        behVals(:,seqLength*2+1) = histcounts(inSeqPresTimes, [0 tsVect])' - histcounts(outSeqPresTimes,[0 tsVect])';
    else
        behVals(:,seqLength*2+1) = behVals(:,seqLength*2+1) + histcounts(inSeqPresTimes, [0 tsVect])' - histcounts(outSeqPresTimes,[0 tsVect])';
    end
end
behDataHeaders{seqLength*2+1} = 'InSeqLog';
itmPresTimes = [ssnData.ItemPresentationTime];
trialPerformance = [ssnData.Performance];
corrTrials = itmPresTimes(logical(trialPerformance));
corTrlHistCounts = histcounts(corrTrials, [0 tsVect])';
inCorrTrials = itmPresTimes(~logical(trialPerformance));
inCorTrlHistCounts = histcounts(inCorrTrials, [0 tsVect])';
behVals(:,seqLength*2+2) = corTrlHistCounts + (inCorTrlHistCounts*-1);
behDataHeaders{seqLength*2+2} = 'PerformanceLog';

behVals(:,seqLength*2+3) = histcounts([ssnData.OdorTrigPokeTime], [0 tsVect])' - histcounts([ssnData.OdorPokeWithdrawTime], [0 tsVect])';
behDataHeaders{seqLength*2+3} = 'PokeEvents';

behVals(:,seqLength*2+4) = histcounts([ssnData.FrontRewardTime], [0 tsVect])';
behDataHeaders{seqLength*2+4} = 'FrontReward';

behVals(:,seqLength*2+5) = histcounts([ssnData.BackRewardTime], [0 tsVect])';
behDataHeaders{seqLength*2+5} = 'BackReward';

if isfield(behaviorData.Raw, 'ErrorSignalTime')
    behVals(:,seqLength*2+6) = histcounts([ssnData.ErrorSignalTime], [0 tsVect])';
    behDataHeaders{seqLength*2+6} = 'ErrorSignal';
end

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
            if isfield(behaviorData.Raw, 'ErrorSignalTime')
                behVals(aniPosHistBins,seqLength*2+7) = aniPosition(:,2);
                behVals(aniPosHistBins,seqLength*2+8) = aniPosition(:,3);
            else
                behVals(aniPosHistBins,seqLength*2+6) = aniPosition(:,2);
                behVals(aniPosHistBins,seqLength*2+7) = aniPosition(:,3);
            end
            findingStrobed = 0;
        end
    end
end
if isfield(behaviorData.Raw, 'ErrorSignalTime')
    behDataHeaders{seqLength*2+7} = 'XvalRatMazePosition';
    behDataHeaders{seqLength*2+8} = 'YvalRatMazePosition';
else
    behDataHeaders{seqLength*2+6} = 'XvalRatMazePosition';
    behDataHeaders{seqLength*2+7} = 'YvalRatMazePosition';
end
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
    tetUnitIDs = cell(1,curNumUnis);
    meanTemplates = repmat({cell(1,4)}, [1,curNumUnis]);
    meanTemplateValues = repmat({nan(3,2)}, [1, curNumUnis]);
    stdTemplates = repmat({cell(1,4)}, [1,curNumUnis]);
    meanSpkRt = cell(1,curNumUnis);
    spkWdth = cell(1,curNumUnis);
    phaseVals = struct('Mean', nan, 'Median', nan, 'R_Length', nan,...
        'AngDev', nan, 'CircStDev', nan, 'R_Test', [nan, nan]);
    phaseInfo = repmat({struct('Theta', phaseVals, 'LowBeta', phaseVals,...
        'Beta', phaseVals, 'LowGamma', phaseVals, 'HighGamma', phaseVals,...
        'Ripple', phaseVals)}, [1,curNumUnis]);
    for uni = 1:curNumUnis
        [~, tempTSs] = plx_ts(curSpkFile, 1, uni);
        tetUnitTSs{uni} = histcounts(tempTSs, [0 tsVect])';
        tetUnitIDs{uni} = [curTet '-U' num2str(uni)];
        for chan = 1:4
            [~, ~, ~, wave] = plx_waves_v(curSpkFile, chan, uni);
            meanTemplates{uni}{chan} = mean(wave);
            stdTemplates{uni}{chan} = std(wave);
        end
        % Calculate the spike width using the spike template
        % First find the largest template of the four wires
        templateHeights = cellfun(@(a)diff([min(a) max(a)]), meanTemplates{uni});
        largestTemplateNum = find(max(templateHeights)==templateHeights,1,'first');
        largestTemplate = meanTemplates{uni}{largestTemplateNum};
        meanTemplateValues{uni}(1,1) = largestTemplateNum;
        % Find the valley
        valleyNdx = find(largestTemplate==min(largestTemplate),1,'first');
        meanTemplateValues{uni}(2,:) = [valleyNdx, largestTemplate(valleyNdx)];
        % Bag the peak... following the valley
        peakNdx = find(largestTemplate(valleyNdx:end)==max(largestTemplate(valleyNdx:end)),1,'first') + (valleyNdx-1);
        meanTemplateValues{uni}(3,:) = [peakNdx, largestTemplate(peakNdx)];
        % Determine the width
        spkWdth{uni} = (peakNdx-valleyNdx)/samp;
        
        % Calculate the overall spike rate (spk/s)
        meanSpkRt{uni} = length(tempTSs)/(size(tetUnitTSs{uni},1)/samp);
        
        % Now pull out summary phase values
        % Theta
        thetaUniPhase = tetThetaHilb(tetUnitTSs{uni}==1);
        phaseInfo{uni}.Theta.Mean = circ_mean(thetaUniPhase);
        phaseInfo{uni}.Theta.Median = circ_median(thetaUniPhase);
        phaseInfo{uni}.Theta.R_Length = circ_r(thetaUniPhase);
        [phaseInfo{uni}.Theta.AngDev, phaseInfo{uni}.Theta.CircStDev] = circ_std(thetaUniPhase);
        [phaseInfo{uni}.Theta.R_Test(1), phaseInfo{uni}.Theta.R_Test(2)] = circ_rtest(thetaUniPhase);
        % Low Beta
        lowBetaUniPhase = tetLowBetaHilb(tetUnitTSs{uni}==1);
        phaseInfo{uni}.LowBeta.Mean = circ_mean(lowBetaUniPhase);
        phaseInfo{uni}.LowBeta.Median = circ_median(lowBetaUniPhase);
        phaseInfo{uni}.LowBeta.R_Length = circ_r(lowBetaUniPhase);
        [phaseInfo{uni}.LowBeta.AngDev, phaseInfo{uni}.LowBeta.CircStDev] = circ_std(lowBetaUniPhase);
        [phaseInfo{uni}.LowBeta.R_Test(1), phaseInfo{uni}.LowBeta.R_Test(2)] = circ_rtest(lowBetaUniPhase);
        % Beta
        betaUniPhase = tetBetaHilb(tetUnitTSs{uni}==1);
        phaseInfo{uni}.Beta.Mean = circ_mean(betaUniPhase);
        phaseInfo{uni}.Beta.Median = circ_median(betaUniPhase);
        phaseInfo{uni}.Beta.R_Length = circ_r(betaUniPhase);
        [phaseInfo{uni}.Beta.AngDev, phaseInfo{uni}.Beta.CircStDev] = circ_std(betaUniPhase);
        [phaseInfo{uni}.Beta.R_Test(1), phaseInfo{uni}.Beta.R_Test(2)] = circ_rtest(betaUniPhase);
        % Low Gamma
        lowGammaUniPhase = tetLowGammaHilb(tetUnitTSs{uni}==1);
        phaseInfo{uni}.LowGamma.Mean = circ_mean(lowGammaUniPhase);
        phaseInfo{uni}.LowGamma.Median = circ_median(lowGammaUniPhase);
        phaseInfo{uni}.LowGamma.R_Length = circ_r(lowGammaUniPhase);
        [phaseInfo{uni}.LowGamma.AngDev, phaseInfo{uni}.LowGamma.CircStDev] = circ_std(lowGammaUniPhase);
        [phaseInfo{uni}.LowGamma.R_Test(1), phaseInfo{uni}.LowGamma.R_Test(2)] = circ_rtest(lowGammaUniPhase);
        % High Gamma
        highGammaUniPhase = tetHighGammaHilb(tetUnitTSs{uni}==1);
        phaseInfo{uni}.HighGamma.Mean = circ_mean(highGammaUniPhase);
        phaseInfo{uni}.HighGamma.Median = circ_median(highGammaUniPhase);
        phaseInfo{uni}.HighGamma.R_Length = circ_r(highGammaUniPhase);
        [phaseInfo{uni}.HighGamma.AngDev, phaseInfo{uni}.HighGamma.CircStDev] = circ_std(highGammaUniPhase);
        [phaseInfo{uni}.HighGamma.R_Test(1), phaseInfo{uni}.HighGamma.R_Test(2)] = circ_rtest(highGammaUniPhase);
        % Ripple
        rippleUniPhase = tetRipHilb(tetUnitTSs{uni}==1);
        phaseInfo{uni}.Ripple.Mean = circ_mean(rippleUniPhase);
        phaseInfo{uni}.Ripple.Median = circ_median(rippleUniPhase);
        phaseInfo{uni}.Ripple.R_Length = circ_r(rippleUniPhase);
        [phaseInfo{uni}.Ripple.AngDev, phaseInfo{uni}.Ripple.CircStDev] = circ_std(rippleUniPhase);
        [phaseInfo{uni}.Ripple.R_Test(1), phaseInfo{uni}.Ripple.R_Test(2)] = circ_rtest(rippleUniPhase);
    end
    unitSummary = struct('UnitName', tetUnitIDs, 'Spike_Features', meanTemplateValues,...
        'TemplateMean', meanTemplates, 'TemplateStDev', stdTemplates,...
        'Mean_SpikeRate', meanSpkRt, 'Spike_Width', spkWdth,...
        'Spike_Phase_Relations', phaseInfo); %#ok<NASGU>
    statMatrix = [tsVect',tetRawLFP, tetRawHilb,...
        tetThetaLFP, tetThetaHilb,...
        tetLowBetaLFP, tetLowBetaHilb,...
        tetBetaLFP, tetBetaHilb,...          
        tetLowGammaLFP, tetLowGammaHilb,...
        tetHighGammaLFP, tetHighGammaHilb,...
        tetRipLFP, tetRipHilb,...
        cell2mat(tetUnitTSs), behVals]; %#ok<NASGU>
    statMatrixColIDs = [{'TimeBin'}, {[curTet '_LFP_Raw']}, {[curTet '_LFP_Raw_HilbVals']},...
        {[curTet '_LFP_Theta']}, {[curTet '_LFP_Theta_HilbVals']},...
        {[curTet '_LFP_LowBeta']}, {[curTet '_LFP_LowBeta_HilbVals']},...
        {[curTet '_LFP_Beta']}, {[curTet '_LFP_Beta_HilbVals']},...
        {[curTet '_LFP_LowGamma']}, {[curTet '_LFP_LowGamma_HilbVals']},...
        {[curTet '_LFP_HighGamma']}, {[curTet '_LFP_HighGamma_HilbVals']},...
        {[curTet '_LFP_Ripple']}, {[curTet '_LFP_Ripple_HilbVals']},...
        tetUnitIDs, behDataHeaders]; %#ok<NASGU>
    save([filePath outputFileName{1} '_' curTet '.mat'], 'statMatrix', 'statMatrixColIDs', 'unitSummary');
    disp([tetNames{tet} ' saved!']);
end

disp('Complete')
cd(origCD);
