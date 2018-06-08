function CreateStatMatrixFromPLX2(plxFile)
%% CreateStatMatrixFromPLX2
% Creates statMatrix data files from a .plx file. This version is made to
% work with the CA1 data files recorded by NJF in Boston. 
% 
% To use this file with other .plx data swap out the behavioral analysis
% function (first line of code in "Behavioral Data" section) with the
% appropriate behavioral analysis file for the rig the data was recorded
% on.
%% Check if filename was fed in, if not identify the file to analyze.
if nargin == 0
    [fileName, filePath] = uigetfile('.plx','Identify .PLX File');
    if fileName == 0
        disp('No file selected');
        return
    end
    plxFile = [filePath fileName];
end
outputFileName = inputdlg('Name Output File', 'Filename', 1, {fileName});

%% Behavioral Data
[behaviorData] = SummarizePLXabbr_BOS(plxFile);
[samp, ~, tetTS, fn, ~] = plx_ad_v(plxFile, 2);
% Create timestamp vector
tsVect = ((0:(fn(1)-1))/samp)+tetTS(1);
% Sometimes the file is composed of multiple fragments, in which case these
% should be combined into a single file
for fragNum = 2:length(tetTS)
    tsVect = [tsVect, ((0:(fn(fragNum)-1))/samp)+tetTS(fragNum)]; %#ok<AGROW>
end
seqLength = behaviorData.Summary.SequenceLength;
behVals = nan(length(tsVect), seqLength*2 + 5);
behDataHeaders = cell(1,seqLength*2 + 7);
ssnData = behaviorData.Raw;
% Step through each sequence item/position and identify when odors were
% presented
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


[numChans, chanNames] = plx_event_names(plxFile);
findingStrobed = 1;
while findingStrobed
    for chan = 1:numChans
        curChan = chanNames(chan,:);
        intVals = double(curChan);
        valLim = find(~(intVals==0), 1, 'last');
        strobedChanLog = strcmp(curChan(1:valLim), 'Strobed');
        if strobedChanLog
            [~, strobedTS, strobedSV] = plx_event_ts(plxFile, curChan);
            [~, ~, ~, aniPosition] = plx_vt_interpret(strobedTS, strobedSV);
            % Remove extra position points that aren't present in LFP
            % signal... not sure why this happens... but it does sometimes.
            aniPosition(aniPosition(:,1)>tsVect(end),:) = [];
            aniPosHistBins = find(histcounts(aniPosition(:,1), [0 tsVect]));
            behVals(aniPosHistBins,seqLength+seqLength+6) = aniPosition(:,2);
            behVals(aniPosHistBins,seqLength+seqLength+7) = aniPosition(:,3);
            findingStrobed = 0;
        end
    end
end
behDataHeaders{seq+seqLength+6} = 'XvalRatMazePosition';
behDataHeaders{seq+seqLength+7} = 'YvalRatMazePosition';
behavMatrix = behVals; %#ok<NASGU>
behavMatrixColIDs = behDataHeaders; %#ok<NASGU>
save([filePath outputFileName{1} '_BehaviorMatrix.mat'], 'behavMatrix', 'behavMatrixColIDs');
disp('Behavior data saved.');

%% Neural Data
[tsCountFl, ~, ~, contCountFl] = plx_info(plxFile, 1);
[~, plxADchanNames] = plx_adchan_names(plxFile);
tetChans = 1:4:(size(tsCountFl,2)-1);
for chan = 2:4:size(tsCountFl,2)
    if sum((chan-1)==tetChans)==1
        curNumUnis = sum(~tsCountFl(2:end, chan)==0);
        tetNum = find(tetChans==(chan-1));
        tetName = ['T' num2str(tetNum)];
        r = 0;
        adChan = [];
        while r < size(plxADchanNames,1)
            r = r + 1;
            chanCheck1a = ~isempty(strfind(plxADchanNames(r,:), [tetName '-1']));
            chanCheck1b = ~isempty(strfind(plxADchanNames(r,:), [tetName '_1']));
            chanCheck2 = ~isempty(strfind(plxADchanNames(r,:), [tetName '-3']));
            chanCheck2b = ~isempty(strfind(plxADchanNames(r,:), [tetName '_3']));
            if contCountFl(r)>0
                if chanCheck1a
                    adChan = [tetName '-1'];
                    break
                elseif chanCheck1b
                    adChan = [tetName '_1'];
                    break
                elseif chanCheck2
                    adChan = [tetName '-3'];
                    break
                elseif chanCheck2b
                    adChan = [tetName '_3'];
                    break
                end
            end
        end
        
        if ~isempty(adChan)
            [samp, ~, tetTS, tetFN, tetV] = plx_ad_v(plxFile, adChan);
            tsVect = ((0:(tetFN(1)-1))/samp)+tetTS(1);
            for fragNum = 2:length(tetTS)
                tsVect = [tsVect, ((0:(tetFN(fragNum)-1))/samp)+tetTS(fragNum)]; %#ok<AGROW>
            end
            [tetRawHilb, tetRawLFP] = PhaseFreqDetectAbbr(tetV,tsVect);
            [tetThetaHilb, tetThetaLFP] = PhaseFreqDetectAbbr(tetV,tsVect,4,12);
            [tetLowBetaHilb, tetLowBetaLFP] = PhaseFreqDetectAbbr(tetV,tsVect,13,19);
            [tetBetaHilb, tetBetaLFP] = PhaseFreqDetectAbbr(tetV,tsVect,20,40);
            [tetLowGammaHilb, tetLowGammaLFP] = PhaseFreqDetectAbbr(tetV,tsVect,41,59);
            [tetHighGammaHilb, tetHighGammaLFP] = PhaseFreqDetectAbbr(tetV,tsVect,61,80);
            [tetRipHilb, tetRipLFP] = PhaseFreqDetectAbbr(tetV,tsVect,150,250);
            tetUnitTSs = cell(1,curNumUnis);
            tetUnitIDs = cell(1,curNumUnis);
            for uni = 1:curNumUnis
                [~, tempTSs] = plx_ts(plxFile, chan-1, uni);
                tetUnitTSs{uni} = histcounts(tempTSs, [0 tsVect])';
                tetUnitIDs{uni} = [tetName '-U' num2str(uni)];
            end
            statMatrix = [tsVect',tetRawLFP, tetRawHilb,...
                tetThetaLFP, tetThetaHilb,...
                tetLowBetaLFP, tetLowBetaHilb,...
                tetBetaLFP, tetBetaHilb,...
                tetLowGammaLFP, tetLowGammaHilb,...
                tetHighGammaLFP, tetHighGammaHilb,...
                tetRipLFP, tetRipHilb,...
                cell2mat(tetUnitTSs), behVals]; %#ok<NASGU>
            statMatrixColIDs = [{'TimeBin'}, {[tetName '_LFP_Raw']}, {[tetName '_LFP_Raw_HilbVals']},...
                {[tetName '_LFP_Theta']}, {[tetName '_LFP_Theta_HilbVals']},...
                {[tetName '_LFP_LowBeta']}, {[tetName '_LFP_LowBeta_HilbVals']},...
                {[tetName '_LFP_Beta']}, {[tetName '_LFP_Beta_HilbVals']},...
                {[tetName '_LFP_LowGamma']}, {[tetName '_LFP_LowGamma_HilbVals']},...
                {[tetName '_LFP_HighGamma']}, {[tetName '_LFP_HighGamma_HilbVals']},...
                {[tetName '_LFP_Ripple']}, {[tetName '_LFP_Ripple_HilbVals']},...
                tetUnitIDs, behDataHeaders]; %#ok<NASGU>
            save([filePath outputFileName{1} '_' tetName '.mat'], 'statMatrix', 'statMatrixColIDs');
            disp([tetName ' saved!']);
        else
            disp(['No ADchan identified for ' tetName]);
        end
    end
end

plx_close(plxFile);
disp('Complete')
