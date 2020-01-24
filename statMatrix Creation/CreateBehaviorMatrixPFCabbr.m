function [behavMatrix, behavMatrixColIDs] = CreateBehaviorMatrixPLXabbr
[plxData] = SummarizePLXevents_SD;
file = plxData.PLXfile;

chan = 1;
[samp, ~, tetTS, fn, ~] = plx_ad_v(file, chan);
while tetTS == -1
    chan = chan + 1;
    [samp, ~, tetTS, fn, ~] = plx_ad_v(file, chan);
end

% Create timestamp vector
tsVect = ((0:(fn(1)-1))/samp)+tetTS(1);
% Sometimes the file is composed of multiple fragments, in which case these
% should be combined into a single file
for fragNum = 2:length(tetTS)
    tsVect = [tsVect, ((0:(fn(fragNum)-1))/samp)+tetTS(fragNum)]; %#ok<AGROW>
end
tsVect(end+1) = tsVect(end)+(1/samp);

seqLength = summary.SequenceLength;
maxSeqLength = max([behaviorData.OrdinalPosition]);
behPad = seqLength + maxSeqLength;
if isfield(summary, 'DualListLog') && summary.DualListLog
    behPad = behPad + seqLength;
end
behVals = nan(length(tsVect)-1, behPad + 5);
behDataHeaders = cell(1,behPad + 7);
% Step through each sequence item/position and identify when odors were
% presented
% Do position first because it doesn't change with multiple lists
% fprintf(outfile, 'Position Counts\n');
for pos = 1:maxSeqLength
    posPresTimes = [behaviorData([behaviorData.OrdinalPosition]==pos).ItemPresentationTime];
    behVals(:,pos) = histcounts(posPresTimes, tsVect)';
    behDataHeaders{pos} = ['Position' num2str(pos)];
%     fprintf(outfile, '     Position #%i = %i trials\n', pos, length(posPresTimes));
end
% fprintf(outfile, 'Odor Counts\n');
for seq = 1:seqLength
    itemPresTimes = [behaviorData([behaviorData.SequenceItem]==seq).ItemPresentationTime];
    behVals(:,seq+maxSeqLength) = histcounts(itemPresTimes, tsVect)';
    behDataHeaders{seq+maxSeqLength} = ['Odor' num2str(seq)];
%     fprintf(outfile, '     Odor #%i = %i trials\n', seq, length(itemPresTimes));
end
if isfield(summary, 'DualListLog') && summary.DualListLog
    for seq = 11:seqLength+10
        itemPresTimes = [behaviorData([behaviorData.SequenceItem]==seq).ItemPresentationTime];
        behVals(:,seq+seqLength+maxSeqLength-10) = histcounts(itemPresTimes, tsVect)';
        behDataHeaders{seq+seqLength+maxSeqLength-10} = ['Odor' num2str(seq)];
%         fprintf(outfile, '     Odor #%i = %i trials\n', seq, length(itemPresTimes));
    end
end

inSeqOdorPres = [behaviorData([behaviorData.TranspositionDistance]==0).ItemPresentationTime];
outSeqOdorPres = [behaviorData(~([behaviorData.TranspositionDistance]==0)).ItemPresentationTime];
behVals(:,behPad+1) = histcounts(inSeqOdorPres, tsVect)' - histcounts(outSeqOdorPres,tsVect)';
behDataHeaders{behPad+1} = 'InSeqLog';
% fprintf(outfile, '\nCompiling InSeq trials.....\n     %i trials were InSeq (%i%%)\n', length(inSeqOdorPres), round(length(inSeqOdorPres)/length(behaviorData),2)*100);
itmPresTimes = [behaviorData.ItemPresentationTime];
trialPerformance = [behaviorData.Performance];
corrTrials = itmPresTimes(logical(trialPerformance));
corTrlHistCounts = histcounts(corrTrials, tsVect)';
inCorrTrials = itmPresTimes(~logical(trialPerformance));
inCorTrlHistCounts = histcounts(inCorrTrials, tsVect)';
behVals(:,behPad+2) = corTrlHistCounts + (inCorTrlHistCounts*-1);
behDataHeaders{behPad+2} = 'PerformanceLog';
% fprintf(outfile, 'Compiling Performance.....\n     %i trials were correct (%i%%)\n', sum(trialPerformance), round(mean(trialPerformance),2)*100);

behVals(:,behPad+3) = histcounts([behaviorData.OdorTrigPokeTime], tsVect)' - histcounts([behaviorData.OdorPokeWithdrawTime], tsVect)';
behDataHeaders{behPad+3} = 'PokeEvents';

behVals(:,behPad+4) = histcounts([behaviorData.FrontRewardTime], tsVect)';
behDataHeaders{behPad+4} = 'FrontReward';

behVals(:,behPad+5) = histcounts([behaviorData.BackRewardTime], tsVect)';
behDataHeaders{behPad+5} = 'BackReward';

if isfield(behaviorData, 'ErrorSignalTime')
    behVals(:,behPad+6) = histcounts([behaviorData.ErrorSignalTime], tsVect)';
    behDataHeaders{behPad+6} = 'ErrorSignal';
end

[numChans, chanNames] = plx_event_names(file);
findingStrobed = 1;
while findingStrobed
    for chan = 1:numChans
        curChan = chanNames(chan,:);
        intVals = double(curChan);
        valLim = find(~(intVals==0), 1, 'last');
        strobedChanLog = strcmp(curChan(1:valLim), 'Strobed');
        if strobedChanLog
            [~, strobedTS, strobedSV] = plx_event_ts(file, curChan);
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
            if isfield(behaviorData, 'ErrorSignalTime')
                behVals(aniPosHistBins,behPad+7) = aniX';
                behVals(aniPosHistBins,behPad+8) = aniY';
            else
                behVals(aniPosHistBins,behPad+6) = aniX;
                behVals(aniPosHistBins,behPad+7) = aniY;
            end
            findingStrobed = 0;
        end
    end
end
if isfield(behaviorData, 'ErrorSignalTime')
    behDataHeaders{behPad+7} = 'XvalRatMazePosition';
    behDataHeaders{behPad+8} = 'YvalRatMazePosition';
else
    behDataHeaders{behPad+6} = 'XvalRatMazePosition';
    behDataHeaders{behPad+7} = 'YvalRatMazePosition';
end
behavMatrix = [tsVect(1:end-1)', behVals];
behavMatrixColIDs = [{'TimeBin'}, behDataHeaders];
save([outputFileName{1} '_BehaviorMatrix.mat'], 'behavMatrix', 'behavMatrixColIDs');
disp('Behavior data saved.');
% fprintf(outfile, 'Behavior Matrix saved as %s_BehaviorMatrix.mat\n', outputFileName{1});
end