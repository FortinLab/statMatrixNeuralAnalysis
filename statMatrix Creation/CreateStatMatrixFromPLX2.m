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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Add in something that enables creation of logical vectors for
%%%%% different periods of the recording... i.e. pre-Session, Session and
%%%%% post-Session recording periods... probably some kind of input system
%%%%% similar to the code to combine ssnData files but operating on the plx
%%%%% files.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Behavioral Data
% [behaviorData] = SummarizePLXabbr_BOS(plxFile);
[behaviorData] = SummarizePLXevents_SD(plxFile);
chan = 1;
[samp, ~, tetTS, fn, ~] = plx_ad_v(plxFile, chan);
while tetTS == -1
    chan = chan + 1;
    [samp, ~, tetTS, fn, ~] = plx_ad_v(plxFile, chan);
end

% Create timestamp vector
tsVect = ((0:(fn(1)-1))/samp)+tetTS(1);
% Sometimes the file is composed of multiple fragments, in which case these
% should be combined into a single file
for fragNum = 2:length(tetTS)
    tsVect = [tsVect, ((0:(fn(fragNum)-1))/samp)+tetTS(fragNum)]; %#ok<AGROW>
end
tsVect(end+1) = tsVect(end)+(1/samp);
ssnData = behaviorData.Raw;
seqLength = behaviorData.Summary.SequenceLength;
maxSeqLength = max([ssnData.OrdinalPosition]);
behVals = nan(length(tsVect)-1, seqLength + maxSeqLength + 5);
behDataHeaders = cell(1,seqLength + maxSeqLength + 7);
% Step through each sequence item/position and identify when odors were
% presented
for seq = 1:seqLength
    for pos = 1:maxSeqLength
        itemPresTimes = [ssnData([ssnData.SequenceItem]==seq).ItemPresentationTime];
        posPresTimes = [ssnData([ssnData.OrdinalPosition]==pos).ItemPresentationTime];
        behVals(:,seq) = histcounts(itemPresTimes, tsVect)';
        behVals(:,pos+seqLength) = histcounts(posPresTimes, tsVect)';
        behDataHeaders{seq} = ['Odor' num2str(seq)];
        behDataHeaders{pos+seqLength} = ['Position' num2str(pos)];
    end
end
inSeqOdorPres = [ssnData([ssnData.TranspositionDistance]==0).ItemPresentationTime];
outSeqOdorPres = [ssnData(~([ssnData.TranspositionDistance]==0)).ItemPresentationTime];
behVals(:,seqLength+maxSeqLength+1) = histcounts(inSeqOdorPres, tsVect)' - histcounts(outSeqOdorPres,tsVect)';
behDataHeaders{maxSeqLength+seqLength+1} = 'InSeqLog';
itmPresTimes = [ssnData.ItemPresentationTime];
trialPerformance = [ssnData.Performance];
corrTrials = itmPresTimes(logical(trialPerformance));
corTrlHistCounts = histcounts(corrTrials, tsVect)';
inCorrTrials = itmPresTimes(~logical(trialPerformance));
inCorTrlHistCounts = histcounts(inCorrTrials, tsVect)';
behVals(:,seqLength + maxSeqLength+2) = corTrlHistCounts + (inCorTrlHistCounts*-1);
behDataHeaders{seqLength + maxSeqLength+2} = 'PerformanceLog';

behVals(:,seqLength + maxSeqLength+3) = histcounts([ssnData.OdorTrigPokeTime], tsVect)' - histcounts([ssnData.OdorPokeWithdrawTime], tsVect)';
behDataHeaders{seqLength + maxSeqLength+3} = 'PokeEvents';

behVals(:,seqLength + maxSeqLength+4) = histcounts([ssnData.FrontRewardTime], tsVect)';
behDataHeaders{seqLength + maxSeqLength+4} = 'FrontReward';

behVals(:,seqLength + maxSeqLength+5) = histcounts([ssnData.BackRewardTime], tsVect)';
behDataHeaders{seqLength + maxSeqLength+5} = 'BackReward';

if isfield(behaviorData.Raw, 'ErrorSignalTime')
    behVals(:,seqLength + maxSeqLength+6) = histcounts([ssnData.ErrorSignalTime], tsVect)';
    behDataHeaders{seqLength + maxSeqLength+6} = 'ErrorSignal';
end

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
            aniPosition(aniPosition(:,1)<tsVect(1),:) = [];
            aniX = aniPosition(:,2);
            aniY = aniPosition(:,3);
            if size(aniPosition,2)>=5
                for t = 1:size(aniPosition,1)
                    if (aniPosition(t,2)==0 && aniPosition(t,3)==0) &&...
                            (aniPosition(t,4)>0 && aniPosition(t,5)>0)
                        aniX(t) = aniPosition(t,4);
                        aniY(t) = aniPosition(t,5);
                    elseif (aniPosition(t,4)==0 && aniPosition(t,5)==0) &&...
                            (aniPosition(t,2)>0 && aniPosition(t,3)>0)
                        aniX(t) = aniPosition(t,2);
                        aniY(t) = aniPosition(t,3);
                    elseif (aniPosition(t,2)>0 && aniPosition(t,3)>0) &&...
                            (aniPosition(t,3)>0 && aniPosition(t,5)>0)
                        aniX(t) = mean([aniPosition(t,2) aniPosition(t,4)]);
                        aniY(t) = mean([aniPosition(t,3) aniPosition(t,5)]);
                    elseif (aniPosition(t,2)==0 && aniPosition(t,3)==0) &&...
                            (aniPosition(t,3)==0 && aniPosition(t,5)==0)
                        aniX(t) = 0;
                        aniY(t) = 0;
                    end
                end
            end
            aniPosHistBins = find(histcounts(aniPosition(:,1), tsVect));
            if isfield(behaviorData.Raw, 'ErrorSignalTime')
                behVals(aniPosHistBins,seqLength + maxSeqLength+7) = aniX';
                behVals(aniPosHistBins,seqLength + maxSeqLength+8) = aniY';
            else
                behVals(aniPosHistBins,seqLength + maxSeqLength+6) = aniX;
                behVals(aniPosHistBins,seqLength + maxSeqLength+7) = aniY;
            end
            findingStrobed = 0;
        end
    end
end
if isfield(behaviorData.Raw, 'ErrorSignalTime')
    behDataHeaders{seqLength + maxSeqLength+7} = 'XvalRatMazePosition';
    behDataHeaders{seqLength + maxSeqLength+8} = 'YvalRatMazePosition';
else
    behDataHeaders{seqLength + maxSeqLength+6} = 'XvalRatMazePosition';
    behDataHeaders{seqLength + maxSeqLength+7} = 'YvalRatMazePosition';
end
behavMatrix = [tsVect(1:end-1)', behVals]; %#ok<NASGU>
behavMatrixColIDs = [{'TimeBin'}, behDataHeaders]; %#ok<NASGU>
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
        tetName1 = sprintf('T%i', tetNum);
        tetName2 = sprintf('T%02i', tetNum);
        r = 0;
        adChan = [];
        while r < size(plxADchanNames,1)
            r = r + 1;
            chanCheck1a1 = ~isempty(strfind(plxADchanNames(r,:), sprintf('%s-1', tetName1)));
            chanCheck1a2 = ~isempty(strfind(plxADchanNames(r,:), sprintf('%s-1', tetName2)));
            chanCheck1b1 = ~isempty(strfind(plxADchanNames(r,:), sprintf('%s_1', tetName1)));
            chanCheck1b2 = ~isempty(strfind(plxADchanNames(r,:), sprintf('%s_1', tetName2)));
            chanCheck2a1 = ~isempty(strfind(plxADchanNames(r,:), sprintf('%s-3', tetName1)));
            chanCheck2a2 = ~isempty(strfind(plxADchanNames(r,:), sprintf('%s-3', tetName2)));
            chanCheck2b1 = ~isempty(strfind(plxADchanNames(r,:), sprintf('%s_3', tetName1)));
            chanCheck2b2 = ~isempty(strfind(plxADchanNames(r,:), sprintf('%s_3', tetName2)));
            if contCountFl(r)>0
                if chanCheck1a1 || chanCheck1b1 || chanCheck2a1 || chanCheck2b1
                    tetName = tetName1;
                elseif chanCheck1a2 || chanCheck1b2 || chanCheck2a2 || chanCheck2b2
                    tetName = tetName2;
                end
                if chanCheck1a1 || chanCheck1a2
                    adChan = [tetName '-1'];
                    break
                elseif chanCheck1b1 || chanCheck1b2
                    adChan = [tetName '_1'];
                    break
                elseif chanCheck2a1 || chanCheck2a2
                    adChan = [tetName '-3'];
                    break
                elseif chanCheck2b1 || chanCheck2b2
                    adChan = [tetName '_3'];
                    break
                end
            end
        end
        
        if ~isempty(adChan)
            fprintf('Found LFP data for %s on %s...........', tetName, adChan);
            [samp, ~, tetTS, tetFN, tetV] = plx_ad_v(plxFile, adChan);
            tsVect = ((0:(tetFN(1)-1))/samp)+tetTS(1);
            for fragNum = 2:length(tetTS)
                tsVect = [tsVect, ((0:(tetFN(fragNum)-1))/samp)+tetTS(fragNum)]; %#ok<AGROW>
            end
            tsVect(end+1) = tsVect(end)+(1/samp); %#ok<AGROW>

            [tetRawHilb, tetRawLFP] = PhaseFreqDetectAbbr(tetV,tsVect(1:end-1));
            [tetThetaHilb, tetThetaLFP] = PhaseFreqDetectAbbr(tetV,tsVect(1:end-1),4,12);
            [tetLowBetaHilb, tetLowBetaLFP] = PhaseFreqDetectAbbr(tetV,tsVect(1:end-1),13,19);
            [tetBetaHilb, tetBetaLFP] = PhaseFreqDetectAbbr(tetV,tsVect(1:end-1),20,40);
            [tetLowGammaHilb, tetLowGammaLFP] = PhaseFreqDetectAbbr(tetV,tsVect(1:end-1),41,59);
            [tetHighGammaHilb, tetHighGammaLFP] = PhaseFreqDetectAbbr(tetV,tsVect(1:end-1),61,80);
            [tetRipHilb, tetRipLFP] = PhaseFreqDetectAbbr(tetV,tsVect(1:end-1),150,250);
            tetUnitTSs = cell(1,curNumUnis);
            tetUnitIDs = cell(1,curNumUnis);
            for uni = 1:curNumUnis
                [~, tempTSs] = plx_ts(plxFile, chan-1, uni);
                tetUnitTSs{uni} = histcounts(tempTSs, tsVect)';
                tetUnitIDs{uni} = [tetName '-U' num2str(uni)];
            end
            statMatrix = [tsVect(1:end-1)',tetRawLFP, tetRawHilb,...
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
            fprintf('%s saved!\n', tetName);
        else
            fprintf('No AD chan found for %s or %s\n', tetName1, tetName2);
        end
    end
end

plx_close(plxFile);
disp('Complete')
