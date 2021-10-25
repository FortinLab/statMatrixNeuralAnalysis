%% Grab and process the behavioral Data
[plxData] = SummarizePLXevents_SD;
load(plxData.Summary.MATfile);

% Pull together the behavior data into a structure organized by trial
osmBehavior = struct('TrialNum', repmat({nan}, [1,length(ssnData)]),...
    'TrialPosition', repmat({nan}, [1,length(ssnData)]),...
    'TrialOdor', repmat({nan}, [1,length(ssnData)]),...
    'PortEntryTime', repmat({nan}, [1,length(ssnData)]),...
    'OdorDeliveryTime', repmat({nan}, [1,length(ssnData)]),...
    'PortExitTime', repmat({nan}, [1,length(ssnData)]),...
    'PerformanceLog', repmat({false}, [1,length(ssnData)]),...
    'RewardSignalTime', repmat({false}, [1,length(ssnData)]),...
    'ErrorSignalTime', repmat({false}, [1,length(ssnData)]),...
    'QuestionableTrialLog', repmat({false}, [1,length(ssnData)]),...
    'BehaviorPLX', cell(1,length(ssnData)),...
    'BehaviorMAT', cell(1,length(ssnData)));

for trl = 1:length(ssnData)
    % Pull out and verify Behavior Data
    curPLX = plxData.Raw(trl);
    curMAT = ssnData(trl);
    osmBehavior(trl).TrialNum = trl;
    if curPLX.OrdinalPosition == curMAT.TrialPosition
        osmBehavior(trl).TrialPosition = curPLX.OrdinalPosition;
    else
        error('Trial position values do not match');
    end
    if curPLX.SequenceItem == curMAT.Odor
        osmBehavior(trl).TrialOdor = curPLX.SequenceItem;
    else
        error('Trial odor values do not match');
    end
    if curPLX.Performance == curMAT.Performance
        osmBehavior(trl).PerformanceLog = curPLX.Performance;
    else
        error('Trial performance values do not match');
    end
    osmBehavior(trl).PortEntryTime = curPLX.OdorTrigPokeTime;
    osmBehavior(trl).OdorDeliveryTime = curPLX.ItemPresentationTime;
    osmBehavior(trl).PortExitTime = curPLX.OdorPokeWithdrawTime;
    osmBehavior(trl).RewardSignalTime = curPLX.RewardSignalTime;
    osmBehavior(trl).ErrorSignalTime = curPLX.ErrorSignalTime;
    osmBehavior(trl).FrontRewardTime = curPLX.FrontRewardTime;
    osmBehavior(trl).BackRewardTime = curPLX.BackRewardTime;
    osmBehavior(trl).QuestionableTrialLog = curPLX.QuestionableTrialLog;
    osmBehavior(trl).BehaviorPLX = curPLX;
    osmBehavior(trl).BehaviorMAT = curMAT;
end

%% Grab Neural data
% First pull out the LFP channel IDs
[~, plxADchanNames] = plx_adchan_names(plxData.Summary.PLXfile);
tetLFPchanNames = cell(size(plxADchanNames,1),1);
for tet = 1:size(plxADchanNames,1)
    tetLFPchanNames{tet} = deblank(plxADchanNames(tet,:));
end
% Now pull out the spike channel IDs
[~,plxTETchanNames] = plx_chan_names(plxData.Summary.PLXfile);
tetSPKchanNames = cell(size(plxTETchanNames,1),1);
for chan = 1:size(plxTETchanNames,1)
    tetSPKchanNames{chan} = deblank(plxTETchanNames(chan,:));
end
% Now extract the counts of spiking & LFP data then identify tetrode IDs
[~, wfCountFl, ~, contCountFl] = plx_info(plxData.Summary.PLXfile, 1);
lfpDataLog = contCountFl ~= 0;
tetLFPchanNames(~lfpDataLog) = [];
tetLFPchanNames(cellfun(@(a)isempty(a), regexp(tetLFPchanNames, '^T([0-9]*)'))) = [];
tetNames = cellfun(@(a)a(1:end-2), tetLFPchanNames','uniformoutput', 0);

% Pull the LFP and spiking data from each tetrode and put them into a
% common structure variable
osmNeural = struct('TetName', tetNames, 'LFP', cell(1,length(tetNames)),...
    'Spikes', cell(1,length(tetNames)));
for tet = 1:length(tetLFPchanNames)
    [samp, ~, tetTS, fn, osmNeural(tet).LFP] = plx_ad_v(plxData.Summary.PLXfile, tetLFPchanNames{tet});
    tetChan = find(strcmp(tetSPKchanNames, tetLFPchanNames{tet}));
    spkCounts = wfCountFl(:,tetChan+1);
    uni = 0;
    for u = 1:sum(spkCounts>0)
        if uni == 0
            [~, osmNeural(tet).Spikes.Unsorted] = plx_ts(plxData.Summary.PLXfile, tetChan, uni);
        else
            [~, osmNeural(tet).Spikes.(sprintf('U%i', uni))] = plx_ts(plxData.Summary.PLXfile, tetChan, uni);
        end
        uni = uni+1;
    end
end

% Create Timestamp Vector
tsVect = ((0:(fn(1)-1))/samp)+tetTS(1);
% Sometimes the file is composed of multiple fragments, in which case these
% should be combined into a single file
for fragNum = 2:length(tetTS)
    tsVect = [tsVect, ((0:(fn(fragNum)-1))/samp)+tetTS(fragNum)]; %#ok<AGROW>
end

%% Save it all!
save([plxData.Summary.PLXfile(1:end-4), '_OSMdata.mat'], 'osmNeural', 'osmBehavior', 'tsVect');
fprintf('%s saved!\n', [plxData.Summary.PLXfile(1:end-4), '_OSMdata.mat']);