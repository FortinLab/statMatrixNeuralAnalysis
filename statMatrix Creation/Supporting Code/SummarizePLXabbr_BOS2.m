function [plxData] = SummarizePLXabbr_BOS2(plxFile)
%% SummarizePLXabbr_BOS2
% Uses the string sessionDataFile that includes file name and directory.
% This is a version of SummarizePLXabbr.m designed to work with the data
% collected in the Eichenbaum lab that uses different event channels
%
% Version 2 created after it was recognized the poke durations were not
% accurate in version 1.

%#ok<*UNRCH,*NASGU,*AGROW>
%% Load the txtData file
[fn, fp] = uigetfile('.mat', ['Identify the txtDta file associated with' plxFile]);
if fn==0
    disp('Analysis Cancelled');
    plxData = [];
    return
end
load([fp fn]);
txtOdors = [txtDta.Odor]; 
txtPos = [txtDta.Position];
txtPokeDur = [txtDta.PokeDur];
txtPokeIn = [txtDta.PokeIn];
txtPokeOut = [txtDta.PokeOut];
txtPokeNums = [txtDta.NumPokes];
txtTargDur = [txtDta.TargDur];
txtBuffer = [txtDta.Buffer];

identifier = sum(round(clock,0));
outfile = fopen(sprintf('%s_txtDta_MergeSummary_%i.txt', plxFile(1:end-4), identifier), 'wt');
fprintf(outfile, '%s txtDta file merger summary\n', plxFile);


%% Turn PLX file into structure
% Compile Plexon Channels into a single structure
% This may be a but unnecessary given the hard coding I do immediately
% after this but a) it's easier to work with this structure and b)
% automating the channel data this way allows for a subfunction to work
% with new channel definitions should the event channels change...
% which they probably will when the rig is upgraded to have 10 odors
[numChans, chanNames] = plx_event_names(plxFile);
plxStruct = struct('Channel', [], 'n', [], 'ts', [], 'sv', []);
plxSsn = cell(numChans,1);
for chan = 1:numChans
    curChan = chanNames(chan,:);
    intVals = double(curChan);
    valLim = find(~(intVals==0), 1, 'last');
    fprintf('Extracting %s timestamps (EventChannel %i of %i)\n', curChan(1:valLim), chan, numChans);
    plxStruct(chan).Channel = curChan(1:valLim);
    [plxStruct(chan).n, plxStruct(chan).ts, plxStruct(chan).sv] = plx_event_ts(plxFile, curChan);
    plxSsn{chan} = [plxStruct(chan).ts, ones(size(plxStruct(chan).ts)).*chan];
end
plx_close(plxFile);
channels = {plxStruct.Channel};
strobeVal = {plxStruct.sv};

plxSummary.Errors = [];
plxSummary.TargetDuration = unique([txtDta.TargDur]);
plxSummary.PlxFile = plxFile;
plxSummary.Identifier = identifier;

%% Identify when Sequence Blocks Start
seqStartChannelLog = strcmp(channels, 'StartArea');
sequenceBlockInitiationTimes = plxStruct(seqStartChannelLog).ts;

%% Identify Poke Times
% The key information used here is when the poke was initiated and how
% long it was held. When the poke was broken is less important and
% easily derived from (initiation time) + (hold time).
pokeInChanNum = strcmp(channels, 'Photocell (IN)');
pokeInitiationTimes = plxStruct(pokeInChanNum).ts;
pokeOutChanNum = strcmp(channels, 'Photocell (OUT)');
pokeEndTimes = plxStruct(pokeOutChanNum).ts;
if ~(length(pokeInitiationTimes) == length(pokeEndTimes))
    pokeInNums = [pokeInitiationTimes ones(length(pokeInitiationTimes),1)];
    pokeOutNums = [pokeEndTimes (ones(length(pokeEndTimes),1).*2)];
    allPokes = sortrows([pokeInNums; pokeOutNums]);
    pokeChanDiff = diff(allPokes(:,2));
    pokeChanRep = find(pokeChanDiff==0);
    if isempty(pokeChanRep) && (length(pokeInitiationTimes) > length(pokeEndTimes))
        pokeInitiationTimes(end) = [];
    elseif isempty(pokeChanRep) && (length(pokeInitiationTimes) < length(pokeEndTimes))
        pokeEndTimes(1) = [];
    elseif ~isempty(pokeChanRep)
        if allPokes(pokeChanRep,2)==2
            pokeOutChanNumLog = pokeEndTimes==allPokes(pokeChanRep+1,1);
            pokeEndTimes(pokeOutChanNumLog) = [];
        elseif allPokes(pokeChanRep,2)==1
            pokeInitiationTimeLog = pokeInitiationTimes==allPokes(pokeChanRep,1);
            pokeInitiationTimes(pokeInitiationTimeLog) = [];
        end
    end
elseif length(pokeInitiationTimes)-length(pokeEndTimes) > 1 || length(pokeInitiationTimes)-length(pokeEndTimes) < 0
    error('Poke In and Out do not add up, something is wrong!')
end
pokeDuration = pokeEndTimes - pokeInitiationTimes;

interPokeInterval = diff([[0; pokeEndTimes], [pokeInitiationTimes; 0]],1,2);
interPokeInterval(1) = []; % Remove first inter-poke interval because it should refer to period between first and second poke, not period from start to end of first poke
interPokeInterval(end) = nan; % NaN out the final entry because there is no poke in after the final poke out
if nansum(interPokeInterval<0)>0
    error('Inter Poke Intervals are not all positive, something is wrong');
end

pokeDiffMtx = pokeEndTimes - pokeInitiationTimes';

%% Identify Trial Odors
% Pull out odor presentation times
odorPresTime = cell(1,10); % ASSUMES 10 ODORS; CHANGE IF THIS ASSUMPTION IS NO LONGER CORRECT

odorAchanNum = strcmp(channels, 'Odor 1');
odorAchanTS = plxStruct(odorAchanNum).ts;
odorAlog = plxStruct(odorAchanNum).sv==0;
odorVlog = plxStruct(odorAchanNum).sv==1;
odorPresTime{1} = odorAchanTS(odorAlog);
odorPresTime{6} = odorAchanTS(odorVlog);

odorBchanNum = strcmp(channels, 'Odor 2');
odorBchanTS = plxStruct(odorBchanNum).ts;
odorBlog = plxStruct(odorBchanNum).sv==0;
odorWlog = plxStruct(odorBchanNum).sv==1;
odorPresTime{2} = odorBchanTS(odorBlog);
odorPresTime{7} = odorBchanTS(odorWlog);

odorCchanNum = strcmp(channels, 'Odor 3');
if ~sum(odorCchanNum)==0
    odorCchanTS = plxStruct(odorCchanNum).ts;
    odorClog = plxStruct(odorCchanNum).sv==0;
    odorXlog = plxStruct(odorCchanNum).sv==1;
    odorPresTime{3} = odorCchanTS(odorClog);
    odorPresTime{8} = odorCchanTS(odorXlog);
end

odorDchanNum = strcmp(channels, 'Odor 4');
if ~sum(odorDchanNum)==0
    odorDchanTS = plxStruct(odorDchanNum).ts;
    odorDlog = plxStruct(odorDchanNum).sv==0;
    odorYlog = plxStruct(odorDchanNum).sv==1;
    odorPresTime{4} = odorDchanTS(odorDlog);
    odorPresTime{9} = odorDchanTS(odorYlog);
end

odorEchanNum = strcmp(channels, 'Odor 5');
if ~sum(odorEchanNum)==0
    odorEchanTS = plxStruct(odorEchanNum).ts;
    odorElog = plxStruct(odorEchanNum).sv==0;
    odorZlog = plxStruct(odorEchanNum).sv==1;
    odorPresTime{5} = odorEchanTS(odorElog);
    odorPresTime{10} = odorEchanTS(odorZlog);
end

% Clean up Odor Presentation Times & Determine Odor Delivery and Poke
% Initiation times
odorTimesAll = cellfun(@(a,b)[a(:), ones(length(a),1).*b], odorPresTime(1:5), num2cell(1:5), 'uniformoutput',0);
odorTimesAllSorted = sortrows(cell2mat(odorTimesAll(:)));
culledTxtTrls = 1;
sequenceLength = max(odorTimesAllSorted(:,2));
o = 1;
while o < size(odorTimesAllSorted,1)
    for o = sequenceLength:size(odorTimesAllSorted,1)
        if sum(odorTimesAllSorted(o-(sequenceLength-1):o,2)' == 1:sequenceLength) == sequenceLength &&...
                sum(diff(odorTimesAllSorted(o-(sequenceLength-1):o,1))<1) == (sequenceLength-1)
            odorTimesAllSorted(o-(sequenceLength-1):o,:) = [];
            if isempty(plxSummary.Errors)
                plxSummary.Errors = {sprintf('Captured odor release on trial %i', o-(sequenceLength-1))};
                fprintf(outfile, 'Captured odor release on trial %i\n', o-(sequenceLength-1));
            else
                plxSummary.Errors{end+1} = sprintf('Captured odor release on trial %i', o-(sequenceLength-1));
                fprintf(outfile, 'Captured odor release on trial %i\n', o-(sequenceLength-1));
            end
            break
        end
    end
end

% Remove spurious odor A flags
dblOdrA = find(diff(odorTimesAllSorted(:,1))<0.001);
if sum(odorTimesAllSorted(dblOdrA,2)==1) == length(dblOdrA)
    odorTimesAllSorted(dblOdrA,:) = [];
else
    error('Not all spurious odors are odor A, check this out');
end

while size(odorTimesAllSorted,1) ~= length(txtDta)
    fprintf(outfile, '\n***** Plexon and Trial Odor Order Don''t Match *****\n');
    fprintf(outfile, 'Trial#:   Plexon Odor     Text Odor\n');
    
    for o = 1:max([size(odorTimesAllSorted,1), length(txtDta)])
        if o>size(odorTimesAllSorted,1)
            plxOdor = nan;
        else
            plxOdor = odorTimesAllSorted(o,2);
        end
        if o>length(txtDta)
            txtOdor = nan;
        else
            txtOdor = txtDta(o).Odor;
        end
        fprintf(outfile, '   %03i          %i          %i\n', o, plxOdor, txtOdor);
    end
    error('Text and Plexon odor order don''t match, clean up the txtDta file.')
end

%% Identify reward times
% To make things easier I'm using the first activation of the solenoid
% as the reward presentation time. Since the solenoid was activated 4x
% times every time reward was presented it's easy to identify the
% initial activation of the solenoid by looking at the inter-solenoid
% interval. Since trials were separated by at least 1.8s any solenoid
% activity with a inter-solenoid interval of <1s has to be related to
% reward presentation.
%
% This is done for front and back reward.

% Front Reward
waterAchanNum = strcmp(channels, 'Water (front)');
interSolenoidInterval = diff([0; plxStruct(waterAchanNum).ts]); % Identify how long it was in between solenoid activations. The zero here is just so that in the following steps 1 is the first spot.
firstSolenoidActivationTime = interSolenoidInterval>=1; % The water solenoid is activated 4x times when reward is triggered with <0.5s in between activations. This line finds the times the solenoid when >=1s between activations, i.e. in between reward presentations
frontRewardTimes = plxStruct(waterAchanNum).ts(firstSolenoidActivationTime);
clear interSolenoidInterval firstSolenoidActivationTime

% Rear Reward
waterBchanNum = strcmp(channels, 'Water (back)');
interSolenoidInterval = diff([0; plxStruct(waterBchanNum).ts]); % Identify how long it was in between solenoid activations. The zero here is just so that in the following steps 1 is the first spot.
firstSolenoidActivationTime = interSolenoidInterval>=1; % The water solenoid is activated 4x times when reward is triggered with <0.5s in between activations. This line finds the times the solenoid when >=1s between activations, i.e. in between reward presentations
backRewardTimes = plxStruct(waterBchanNum).ts(firstSolenoidActivationTime);
clear interSolenoidInterval firstSolenoidActivationTime
%
% if sum([txtDta.Perf]==1 & [txtDta.TransDist]==0) == length(frontRewardTimes)
rwdCase = 'InSeq';
% elseif sum([txtDta.Perf]==1) == length(frontRewardTimes)
%     rwdCase = 'Both';
% else
%     error('Text & Plexon reward times don''t match');
% end

%% Create plxData Structure
% The structure plxData is a 1 x n structure where n = the number of trials
plxSession = struct('TrialNum', {txtDta.TrialNum},...
    'OrdinalPosition', {txtDta.Position}, 'SequenceItem', {txtDta.Odor},...
    'TranspositionDistance', {txtDta.TransDist}, 'SessionBlockNumber', {txtDta.SeqNum},...
    'Performance', {txtDta.Perf},...
    'OdorDeliveryDelay',  repmat({nan}, size(txtDta)),...
    'ItemPresentationTime',  repmat({nan}, size(txtDta)),...
    'OdorTrigPokeTime',  repmat({nan}, size(txtDta)),...
    'PokeDuration', repmat({nan}, size(txtDta)),...
    'OdorPokeWithdrawTime', repmat({nan}, size(txtDta)),...
    'OdorPokesDurations', repmat({[]}, size(txtDta)),...
    'TrialPokesDurations', repmat({[]}, size(txtDta)),...
    'PositionData', repmat({nan}, size(txtDta)),...
    'AniReturnToOrigin', repmat({nan}, size(txtDta)),...
    'FrontRewardTime', repmat({nan}, size(txtDta)),...
    'BackRewardTime', repmat({nan}, size(txtDta)),...
    'MultiOdorPokeLog', repmat({nan}, size(txtDta)),...
    'Errors', repmat({[]}, size(txtDta)));
plxSummary.SequenceLength = sequenceLength;
plxSummary.FileName = plxFile;

%% Run through the session
for trl = 1:length(plxSession)
    clear curTrlPokeMtrx
    curTargDur = txtTargDur(trl);
    curTxtBuffer = txtBuffer(trl)/1000;
    
    plxSession(trl).ItemPresentationTime = odorTimesAllSorted(trl,1);
    plxSession(trl).OdorTrigPokeTime =  pokeInitiationTimes(find(pokeInitiationTimes<odorTimesAllSorted(trl,1), 1, 'last'));
    plxSession(trl).OdorDeliveryDelay = plxSession(trl).ItemPresentationTime - plxSession(trl).OdorTrigPokeTime;
    
    if trl == length(plxSession)
        frontRwdTime = frontRewardTimes(find(frontRewardTimes>odorTimesAllSorted(trl,1),1,'first'));
        rearRwdTime = backRewardTimes(find(backRewardTimes>odorTimesAllSorted(trl,1),1,'first'));
        allTrlPokesLog = pokeInitiationTimes>=plxSession(trl).OdorTrigPokeTime;
    else
        frontRwdTime = frontRewardTimes((frontRewardTimes>odorTimesAllSorted(trl,1)) & (frontRewardTimes<odorTimesAllSorted(trl+1,1)));
        rearRwdTime = backRewardTimes((backRewardTimes>odorTimesAllSorted(trl,1)) & (backRewardTimes<odorTimesAllSorted(trl+1,1)));
        allTrlPokesLog = (pokeInitiationTimes>=plxSession(trl).OdorTrigPokeTime) &...
            (pokeInitiationTimes<odorTimesAllSorted(trl+1,1));
    end
    plxSession(trl).FrontRewardTime = frontRwdTime;
    plxSession(trl).BackRewardTime = rearRwdTime;
        
    allTrlPokesMtrx = [pokeInitiationTimes(allTrlPokesLog),...
        pokeEndTimes(allTrlPokesLog),...
        pokeInitiationTimes(allTrlPokesLog)-plxSession(trl).OdorTrigPokeTime,...
        pokeEndTimes(allTrlPokesLog)-plxSession(trl).OdorTrigPokeTime,...
        interPokeInterval(allTrlPokesLog)];
    if trl ~= length(plxSession)
        allTrlPokesMtrx(end,:) = [];
    end
    plxSession(trl).TrialPokesDurations = allTrlPokesMtrx;
    
%     if (frontRwdTime - odorTimesAllSorted(trl,1)) > (curTargDur + 1)
%         plxSession(trl).Errors = [plxSession(trl).Errors; {sprintf('Excessive reward latency (%i), possible poke in timing error on #%i', round(frontRwdTime - odorTimesAllSorted(trl,1), 4), trl)}];
%     end
    
        
    % Next identify all potential trial pokes
    switch rwdCase
        case 'InSeq'
            if plxSession(trl).TranspositionDistance==0 && plxSession(trl).Performance==1
                if ~isempty(frontRwdTime)
                    potTrlPokesLog = (pokeInitiationTimes>=plxSession(trl).OdorTrigPokeTime) &...
                        (pokeInitiationTimes<frontRwdTime);
                else
                    plxSession(trl).Errors = [plxSession(trl).Errors; {'ErrType1: No reward time stamp when one should be present'}];
                    fprintf(outfile, 'ErrType1: No reward time stamp when one should be present on trial %i\n', trl);
                    potTrlPokesLog = (pokeInitiationTimes>=plxSession(trl).OdorTrigPokeTime) &...
                        (pokeInitiationTimes<(plxSession(trl).ItemPresentationTime + curTargDur));
                end
            else
                potTrlPokesLog = (pokeInitiationTimes>=plxSession(trl).OdorTrigPokeTime) &...
                    (pokeInitiationTimes<(plxSession(trl).ItemPresentationTime + curTargDur));
            end
        case 'Both'
    end
    
    curTrlPokeMtrx = [pokeInitiationTimes(potTrlPokesLog),...
        pokeEndTimes(potTrlPokesLog),...
        pokeInitiationTimes(potTrlPokesLog)-plxSession(trl).OdorTrigPokeTime,...
        pokeEndTimes(potTrlPokesLog)-plxSession(trl).OdorTrigPokeTime,...
        interPokeInterval(potTrlPokesLog)];
    
    plxSession(trl).OdorPokesDurations = curTrlPokeMtrx;
    
    
    if size(curTrlPokeMtrx,1)==1
        plxSession(trl).PokeDuration = curTrlPokeMtrx(1,4);
        plxSession(trl).OdorPokeWithdrawTime = curTrlPokeMtrx(1,2);
        plxSession(trl).MultiOdorPokeLog = false;
        if (round(curTrlPokeMtrx(1,4),4) == round(txtPokeDur(trl),4)) ||...
                round(curTrlPokeMtrx(1,4),5) == round(txtPokeDur(trl),5)
        elseif ((round(curTrlPokeMtrx(1,4),4) > round(txtPokeDur(trl),4))||...
                (round(curTrlPokeMtrx(1,4),5) > round(txtPokeDur(trl),5))) &&...
                ((plxSession(trl).TranspositionDistance==0 && plxSession(trl).Performance==1) ||...
                (plxSession(trl).TranspositionDistance~=0 && plxSession(trl).Performance~=1))
        else
            plxSession(trl).Errors = [plxSession(trl).Errors; {'ErrType2: Here''s a stupid trial'}];
            fprintf(outfile, 'ErrType2: Here''s a stupid trial %i is stupid, probable flick buffer error\n', trl);
        end
    elseif size(curTrlPokeMtrx,1)>=2
        plxSession(trl).MultiOdorPokeLog = true;
        p = 1;
        curIPI = curTrlPokeMtrx(p,5);
        curPokeDur = curTrlPokeMtrx(p,4);
        if curIPI<curTxtBuffer
            while curIPI<curTxtBuffer && p<=size(curTrlPokeMtrx,1)
                curIPI = curTrlPokeMtrx(p,5);
                curPokeDur = curTrlPokeMtrx(p,4);
                if (curIPI > curTxtBuffer) && (p<size(curTrlPokeMtrx,1))
                    plxSession(trl).PokeDuration = curTrlPokeMtrx(end,4);
                    plxSession(trl).OdorPokeWithdrawTime = curTrlPokeMtrx(end,2);
                    plxSession(trl).Errors = [plxSession(trl).Errors; {sprintf('ErrType3: Buffer Error, IPI between pokes %i and %i are longer than buffer', p, p+1)}];
                    fprintf(outfile, 'ErrType3: Buffer Error, IPI between pokes %i and %i are longer than buffer on trial #%i\n', p, p+1, trl);
                    break
                end
                if curPokeDur > curTargDur
                    if (plxSession(trl).TranspositionDistance==0 && plxSession(trl).Performance==1) ||...
                            (plxSession(trl).TranspositionDistance~=0 && plxSession(trl).Performance~=1)
                        plxSession(trl).PokeDuration = curPokeDur;
                        plxSession(trl).OdorPokeWithdrawTime = curTrlPokeMtrx(p,2);
                        break
                    else
                        if txtPokeDur(trl)<curTargDur
                            plxSession(trl).Errors = [plxSession(trl).Errors; {sprintf('ErrType5: Target duration passed when it shouldn''t have on trial #%i', trl)}];
                            fprintf(outfile, '\nErrType5: Target duration passed when it shouldn''t have on trial #%i\n', trl);
                            
                            [out, in] = find((round(txtPokeDur(trl),4)==round(pokeEndTimes - pokeInitiationTimes',4)) |...
                                (round(txtPokeDur(trl),5)==round(pokeEndTimes - pokeInitiationTimes',5)));
                            if ~isempty(in) && pokeInitiationTimes(in)<curTrlPokeMtrx(1,1)
                                error('here');
                            end
                            if ~isempty(out)
                                trueTrlEndTime = pokeEndTimes(out);
                                newPotTrlPokesLog = (pokeInitiationTimes>=plxSession(trl).OdorTrigPokeTime) &...
                                    (pokeInitiationTimes<(trueTrlEndTime+0.01));
                                newCurTrlPokeMtrx = [pokeInitiationTimes(newPotTrlPokesLog),...
                                    pokeEndTimes(newPotTrlPokesLog),...
                                    pokeInitiationTimes(newPotTrlPokesLog)-plxSession(trl).OdorTrigPokeTime,...
                                    pokeEndTimes(newPotTrlPokesLog)-plxSession(trl).OdorTrigPokeTime,...
                                    interPokeInterval(newPotTrlPokesLog)];
                                errPokeInNum = find(newCurTrlPokeMtrx(:,1)==pokeInitiationTimes(in));
                                
                                plxSession(trl).Errors = [plxSession(trl).Errors; {sprintf('ErrType6: Poke in time wrong, trial was timed with poke #%i as poke in on trial #%i', errPokeInNum, trl)}];
                                fprintf(outfile, 'ErrType6: Poke in time wrong, trial was timed with poke #%i as poke in on trial #%i\n', errPokeInNum, trl);
                                plxSession(trl).OdorPokesDurations = newCurTrlPokeMtrx;
                                plxSession(trl).PokeDuration = trueTrlEndTime - plxSession(trl).OdorTrigPokeTime;
                                plxSession(trl).OdorPokeWithdrawTime = curTrlPokeMtrx(p,2);
                                break
                            else
                                plxSession(trl).PokeDuration = curPokeDur;
                                plxSession(trl).OdorPokeWithdrawTime = curTrlPokeMtrx(p,2);
                                break
                            end
                        else 
                            error('Code no should go here... this is a sanity check')
                        end
                    end
                elseif curPokeDur < curTargDur
                    if (round(curPokeDur,4) == round(txtPokeDur(trl),4)) ||...
                            (round(curPokeDur,5) == round(txtPokeDur(trl),5))
                        plxSession(trl).PokeDuration = curPokeDur;
                        plxSession(trl).OdorPokeWithdrawTime = curTrlPokeMtrx(p,2);
                        break
                    elseif (round(txtPokeDur(trl),4) == round(curTrlPokeMtrx(p,2)-curTrlPokeMtrx(p,1),4)) ||...
                            (round(txtPokeDur(trl),5) == round(curTrlPokeMtrx(p,2)-curTrlPokeMtrx(p,1),5))
                        plxSession(trl).PokeDuration = curPokeDur;
                        plxSession(trl).OdorPokeWithdrawTime = curTrlPokeMtrx(p,2);
                        plxSession(trl).Errors = [plxSession(trl).Errors; {sprintf('ErrType7: Poke #%i withdrawal treated as poke initiation error on trial #%i', p, trl)}];
                        fprintf(outfile, 'ErrType7: Poke #%i withdrawal treated as poke initiation error on trial #%i\n', p, trl);
                        break
                    elseif p == size(curTrlPokeMtrx,1)
                        [out, in] = find((round(txtPokeDur(trl),4)==round(pokeDiffMtx,4)) |...
                            (round(txtPokeDur(trl),5)==round(pokeDiffMtx,5)));
                        if ~isempty(in) && pokeInitiationTimes(in)<curTrlPokeMtrx(1,1)
                            error('here');
                        end
                        if ~isempty(out)
                            trueTrlEndTime = pokeEndTimes(out);
                            if trl<length(plxSession)
                                if pokeInitiationTimes(in) < odorTimesAllSorted(trl+1,1)
                                    newPotTrlPokesLog = (pokeInitiationTimes>=plxSession(trl).OdorTrigPokeTime) &...
                                        (pokeInitiationTimes<odorTimesAllSorted(trl+1,1));
                                    newTrlPokeMtrx = [pokeInitiationTimes(newPotTrlPokesLog),...
                                        pokeEndTimes(newPotTrlPokesLog),...
                                        pokeInitiationTimes(newPotTrlPokesLog)-plxSession(trl).OdorTrigPokeTime,...
                                        pokeEndTimes(newPotTrlPokesLog)-plxSession(trl).OdorTrigPokeTime,...
                                        interPokeInterval(newPotTrlPokesLog)];
                                    errPokeInNum = find(newTrlPokeMtrx(:,1)==pokeInitiationTimes(in));
                                    
                                    plxSession(trl).PokeDuration = trueTrlEndTime - plxSession(trl).OdorTrigPokeTime;
                                    plxSession(trl).OdorPokesDurations = newTrlPokeMtrx(1:end-1,:);
                                    plxSession(trl).Errors = [plxSession(trl).Errors; {sprintf('ErrType8: Poke in time wrong, trial was timed with poke #%i as poke in on trial #%i', errPokeInNum,trl)}];
                                    fprintf(outfile, 'ErrType8: Poke in time wrong, trial was timed with poke #%i as poke in on trial #%i\n', errPokeInNum,trl);
                                    break
                                else
                                    plxSession(trl).PokeDuration = curTrlPokeMtrx(end,4);
                                    plxSession(trl).Errors = [plxSession(trl).Errors; {'ErrType9: This trial is stupid'}];
                                    fprintf(outfile, 'ErrType9: Trial #%i is stupid, possible buffer error\n', trl);
                                    break
                                end
                            else
                                newPotTrlPokesLog = (pokeInitiationTimes>=plxSession(trl).OdorTrigPokeTime);
                                newTrlPokeMtrx = [pokeInitiationTimes(newPotTrlPokesLog),...
                                    pokeEndTimes(newPotTrlPokesLog),...
                                    pokeInitiationTimes(newPotTrlPokesLog)-plxSession(trl).OdorTrigPokeTime,...
                                    pokeEndTimes(newPotTrlPokesLog)-plxSession(trl).OdorTrigPokeTime,...
                                    interPokeInterval(newPotTrlPokesLog)];
                                errPokeInNum = find(newTrlPokeMtrx(:,1)==pokeInitiationTimes(in));
                                
                                plxSession(trl).PokeDuration = trueTrlEndTime - plxSession(trl).OdorTrigPokeTime;
                                plxSession(trl).OdorPokesDurations = newTrlPokeMtrx(1:end-1,:);
                                plxSession(trl).Errors = [plxSession(trl).Errors; {sprintf('ErrType10: Poke in time wrong, trial was timed with poke #%i as poke in on trial #%i', errPokeInNum,trl)}];
                                    fprintf(outfile, 'ErrType10: Poke in time wrong, trial was timed with poke #%i as poke in on trial #%i\n', errPokeInNum,trl);
                                break
                            end
                        else
                            plxSession(trl).PokeDuration = curTrlPokeMtrx(end,4);
                            plxSession(trl).OdorPokeWithdrawTime = curTrlPokeMtrx(end,2);
                            plxSession(trl).Errors = [plxSession(trl).Errors; {'ErrType11: No matching poke duration for the textfile hold duration value'}];
                            fprintf(outfile, 'ErrType11: No matching poke duration for the textfile hold duration value for trial #%i\n', trl);
                            break
                        end
                    else
                        p = p + 1;
                    end
                end
            end
        else
            plxSession(trl).PokeDuration = curTrlPokeMtrx(end,4);
            plxSession(trl).OdorPokeWithdrawTime = curTrlPokeMtrx(end,2);
            plxSession(trl).Errors = [plxSession(trl).Errors; {sprintf('ErrType12: Initial poke buffer error on trial #%i', trl)}];
            fprintf(outfile, 'ErrType12: Initial poke buffer error on trial... flicking? #%i\n', trl);
        end
    else
        error('No pokes found...huh?');
    end
    if plxSession(trl).PokeDuration > curTargDur
        if (plxSession(trl).TranspositionDistance == 0) && (plxSession(trl).Performance == 0)
            plxSession(trl).Errors = [plxSession(trl).Errors; {sprintf('ErrType13: InSeq trial held past target, marked incorrect on trial #%i', trl)}];
            fprintf(outfile, 'ErrType13: InSeq trial held past target, marked incorrect on trial #%i\n', trl);
        elseif (plxSession(trl).TranspositionDistance ~= 0) && (plxSession(trl).Performance == 1)
            plxSession(trl).Errors = [plxSession(trl).Errors; {sprintf('ErrType14: OutSeq trial held past target, marked correct on trial #%i', trl)}];
            fprintf(outfile, 'ErrType14: OutSeq trial held past target, marked correct on trial #%i\n', trl);
        end
    elseif plxSession(trl).PokeDuration < curTargDur
        if (plxSession(trl).TranspositionDistance ~= 0) && (plxSession(trl).Performance == 0)
            plxSession(trl).Errors = [plxSession(trl).Errors; {sprintf('ErrType15: OutSeq trial withdrawn, marked incorrect on trial #%i', trl)}];
            fprintf(outfile, 'ErrType15: OutSeq trial withdrawn, marked incorrect on trial #%i\n', trl);
        elseif (plxSession(trl).TranspositionDistance == 0) && (plxSession(trl).Performance == 1)
            plxSession(trl).Errors = [plxSession(trl).Errors; {sprintf('ErrType16: InSeq trial withdrawn, marked correct on trial #%i', trl)}];
            fprintf(outfile, 'ErrType16: InSeq trial withdrawn, marked correct on trial #%i\n', trl);
        end
    end
    if ~isempty(plxSession(trl).Errors)
        fprintf(outfile, 'Poke Initiation for Trial #%i = %f\n', trl, plxSession(trl).OdorTrigPokeTime);
        fprintf(outfile, '\n');
    end
end


%% Create poke diff plot
txtHD = [txtDta.PokeDur];
plxHD = [plxSession.PokeDuration];

h = figure;
inSeqLog = [plxSession.TranspositionDistance]==0;
subplot(2,2,1)
scatter(txtHD, plxHD);
hold on;
scatter(txtHD(inSeqLog), plxHD(inSeqLog), 'r');
set(gca, 'xlim', [0 2], 'ylim', [0 2]);
title(plxFile)
xlabel('Text File Durations');
ylabel('Plexon Durations');
subplot(2,2,3:4)
stem(plxHD-txtHD);
hold on;
stem(find(inSeqLog), plxHD(inSeqLog)-txtHD(inSeqLog), 'r');
set(gca, 'xtick', 0:10:300);
grid on

saveas(h, [plxFile(1:end-4) '_PokeDur_Summary.fig']);
%%
plxData.Raw = plxSession;
plxData.Summary = plxSummary;
plxData.TextData = txtDta;
fclose(outfile);
