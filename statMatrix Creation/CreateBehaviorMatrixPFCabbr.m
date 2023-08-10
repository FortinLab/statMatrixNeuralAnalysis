function CreateBehaviorMatrixPFCabbr
[fileName, path] = uigetfile('.plx','Identify .PLX File');
if fileName == 0
    disp('No file selected, analysis cancelled')
    return
end
plxFile = [path fileName];
[path, fileName] = fileparts(plxFile);
path = [path '\'];
cd(path);
outfile = fopen(sprintf('%s_PLXeventSummary.txt', fileName), 'a+');

[plxData] = SummarizePLXevents_SD(plxFile, [], outfile);
%%
chan = 1;
[samp, ~, tetTS, fn, ~] = plx_ad_v(plxData.Summary.PLXfile, chan);
while tetTS == -1
    chan = chan + 1;
    [samp, ~, tetTS, fn, ~] = plx_ad_v(plxData.Summary.PLXfile, chan);
end

%% Create timestamp vector
tsVect = ((0:(fn(1)-1))/samp)+tetTS(1);
% Sometimes the file is composed of multiple fragments, in which case these
% should be combined into a single file
for fragNum = 2:length(tetTS)
    tsVect = [tsVect, ((0:(fn(fragNum)-1))/samp)+tetTS(fragNum)]; %#ok<AGROW>
end
tsVect(end+1) = tsVect(end)+(1/samp);

seqLength = plxData.Summary.SequenceLength;
maxSeqLength = max([plxData.Raw.OrdinalPosition]);
%%
% Step through each sequence item/position and identify when odors were
% presented
% Do position first because it doesn't change with multiple lists
% fprintf(outfile, 'Position Counts\n');
posVals = nan(length(tsVect)-1, seqLength);
posHeaders = cell(1,seqLength);
for pos = 1:maxSeqLength
    posPresTimes = [plxData.Raw([plxData.Raw.OrdinalPosition]==pos).ItemPresentationTime];
    posVals(:,pos) = histcounts(posPresTimes, tsVect)';
    posHeaders{pos} = ['Position' num2str(pos)];
    fprintf(outfile, '     Position #%i = %i trials\n', pos, length(posPresTimes));
end
% fprintf(outfile, 'Odor Counts\n');
seqVals = nan(length(tsVect)-1, seqLength);
seqHeaders = cell(1,seqLength);
for seq = 1:seqLength
    itemPresTimes = [plxData.Raw([plxData.Raw.SequenceItem]==seq).ItemPresentationTime];
    seqVals(:,seq) = histcounts(itemPresTimes, tsVect)';
    seqHeaders{seq} = ['Odor' num2str(seq)];
    fprintf(outfile, '     Odor #%i = %i trials\n', seq, length(itemPresTimes));
end
if isfield(plxData.Summary, 'DualListLog') && plxData.Summary.DualListLog
    dlSeqVals = nan(length(tsVect)-1,seqLength);
    dlSeqHeaders = cell(1,seqLength);
    for seq = 11:seqLength+10
        itemPresTimes = [plxData.Raw([plxData.Raw.SequenceItem]==seq).ItemPresentationTime];
        dlSeqVals(:,seq) = histcounts(itemPresTimes, tsVect)';
        dlSeqHeaders{seq} = ['Odor' num2str(seq)];
        fprintf(outfile, '     Odor #%i = %i trials\n', seq, length(itemPresTimes));
    end
    seqVals = [seqVals, dlSeqVals];
    seqHeaders = [seqHeaders, dlSeqHeaders];
end
% In Seq Vect
inSeqOdorPres = [plxData.Raw([plxData.Raw.TranspositionDistance]==0).ItemPresentationTime];
outSeqOdorPres = [plxData.Raw(~([plxData.Raw.TranspositionDistance]==0)).ItemPresentationTime];
isVals = histcounts(inSeqOdorPres, tsVect)' - histcounts(outSeqOdorPres,tsVect)';
fprintf(outfile, '\nCompiling InSeq trials.....\n     %i trials were InSeq (%i%%)\n', length(inSeqOdorPres), round(length(inSeqOdorPres)/length(plxData.Raw),2)*100);

itmPresTimes = [plxData.Raw.ItemPresentationTime];
trialPerformance = [plxData.Raw.Performance];
corrTrials = itmPresTimes(logical(trialPerformance));
corTrlHistCounts = histcounts(corrTrials, tsVect)';
inCorrTrials = itmPresTimes(~logical(trialPerformance));
inCorTrlHistCounts = histcounts(inCorrTrials, tsVect)';
perfVals = corTrlHistCounts + (inCorTrlHistCounts*-1);
fprintf(outfile, 'Compiling Performance.....\n     %i trials were correct (%i%%)\n', sum(trialPerformance), round(mean(trialPerformance),2)*100);

lightVals = histcounts([plxData.Raw.TrialLightTime], tsVect)';

pokeVals = histcounts([plxData.Raw.OdorTrigPokeTime], tsVect)' - histcounts([plxData.Raw.OdorPokeWithdrawTime], tsVect)';

rwdVals = histcounts([plxData.Raw.FrontRewardTime], tsVect)';

backRwdVals = histcounts([plxData.Raw.BackRewardTime], tsVect)';

errVals = histcounts([plxData.Raw.ErrorSignalTime], tsVect)';

rwdSigVals = histcounts([plxData.Raw.RewardSignalTime], tsVect)';

[numChans, chanNames] = plx_event_names(plxData.Summary.PLXfile);
findingStrobed = 1;
aniPosVals = nan(length(tsVect)-1, 2);
while findingStrobed
    for chan = 1:numChans
        curChan = chanNames(chan,:);
        intVals = double(curChan);
        valLim = find(~(intVals==0), 1, 'last');
        strobedChanLog = strcmp(curChan(1:valLim), 'Strobed');
        if strobedChanLog
            [~, strobedTS, strobedSV] = plx_event_ts(plxData.Summary.PLXfile, curChan);
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
            aniPosVals(aniPosHistBins,1) = aniX';
            aniPosVals(aniPosHistBins,2) = aniY';
            findingStrobed = 0;
        end
    end
end

behavMatrix = [tsVect(1:end-1)', posVals, seqVals, isVals, perfVals, lightVals, pokeVals,...
    rwdVals, backRwdVals, errVals, rwdSigVals, aniPosVals];
behavMatrixColIDs = ['TimeBin', posHeaders, seqHeaders, 'InSeqLog', 'PerformanceLog', 'PortLight', 'PokeEvents',...
    'FrontReward', 'BackReward', 'ErrorSignal', 'RewardSignal', 'XvalRatMazePosition', 'YvalRatMazePosition'];

save([plxData.Summary.PLXfile(1:end-4) '_BehaviorMatrix.mat'], 'behavMatrix', 'behavMatrixColIDs', 'plxData');
disp('Behavior data saved.');
fprintf(outfile, 'Behavior Matrix saved as %s_BehaviorMatrix.mat\n', plxData.Summary.PLXfile(1:end-4));
end