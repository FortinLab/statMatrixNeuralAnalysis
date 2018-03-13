function [plxData] = SummarizePLXabbr_BOS(plxFile)
%% SummarizePLXabbr_BOS
% Uses the string sessionDataFile that includes file name and directory.
% This is a version of SummarizePLXabbr.m designed to work with the data
% collected in the Eichenbaum lab that uses different event channels

%#ok<*UNRCH,*NASGU,*AGROW>
%% Create plxData Structure
% The structure plxData is a 1 x n structure where n = the number of trials
plxSession = struct('OrdinalPosition', [], 'SequenceItem', [], 'TranspositionDistance', [],...
    'SessionBlockNumber', [], 'Performance', [], 'SessionBlockStartTime', [],...
    'ItemPresentationTime', [], 'OdorTrigPokeTime', [], 'PokeDuration', [],...
    'OdorPokeWithdrawTime', [], 'OdorPokesDurations', [], 'PositionData', [], 'AniReturnToOrigin', [],...
    'FrontRewardTime', [], 'BackRewardTime', [], 'SpuriousOdorA', nan,...
    'MultiOdorPokeLog', []);
plxSummary.Errors = [];
%% Turn PLX file into structure
% Compile Plexon Channels into a single structure
% This may be a but unnecessary given the hard coding I do immediately
% after this but a) it's easier to work with this structure and b)
% automating the channel data this way allows for a subfunction to work
% with new channel definitions should the event channels change...
% which they probably will when the rig is upgraded to have 10 odors
[numChans, chanNames] = plx_event_names(plxFile);
plxStruct = struct('Channel', [], 'n', [], 'ts', [], 'sv', []);
for chan = 1:numChans
    curChan = chanNames(chan,:);
    intVals = double(curChan);
    valLim = find(~(intVals==0), 1, 'last');
    plxStruct(chan).Channel = curChan(1:valLim);
    [plxStruct(chan).n, plxStruct(chan).ts, plxStruct(chan).sv] = plx_event_ts(plxFile, curChan);
end
plx_close(plxFile);
channels = {plxStruct.Channel};

%% Pull out "Start" Times
startChannelLog = strcmp(channels, 'Start');
startTimes = plxStruct(startChannelLog).ts;
if isempty(startTimes)
    error('No Start Time!');
elseif length(startTimes) > 1
    multStart = 1;
else
    multStart = 0;
end

%% Check for Activity on the Terminate Channel
% Read the terminate channel and determine if a terminate command was
% issued anytime during the session. If so then truncate all the other
% behavioral variables to only the time points that occurred prior to
% when the command was issued
terminateChannelLog = strcmp(channels, 'Quit Matlab (Alt+9)');
if ~(sum(terminateChannelLog)==0) && ~(plxStruct(terminateChannelLog).n == 0)
    term = 1;
    plxSummary.Terminate = 1;
    terminateTime = plxStruct(terminateChannelLog).ts(find(plxStruct(terminateChannelLog).ts>0, 1, 'first'));
    if multStart
        termTimes = plxStruct(terminateChannelLog).ts;
        behPeriods = [startTimes, nan(size(startTimes))];
        if ~isempty(termTimes(termTimes<startTimes(end) & termTimes>startTimes(1)))
            for start = 2:length(startTimes)
                intermedTermFlags = termTimes(termTimes<startTimes(start) & termTimes>startTimes(start-1));
                behPeriods(start-1,2) = intermedTermFlags(1);
            end
            if sum(termTimes>startTimes(end))>=1
                behPeriods(end,2) = termTimes(find(termTimes>startTimes(end),1,'first'));
            else
                behPeriods(end,2) = max(cell2mat({plxStruct.ts}')) + 10;
            end
        else
            multStart = 0;
        end
    end
else
    term = 0;
    plxSummary.Terminate = 0;
end

%% Identify when Sequence Blocks Start
seqStartChannelLog = strcmp(channels, 'StartArea');
sequenceBlockInitiationTimes = plxStruct(seqStartChannelLog).ts;

if term && ~multStart
    sequenceBlockInitiationTimes(sequenceBlockInitiationTimes>terminateTime) = [];
elseif term && multStart
    seqBlockInitTimesBehPrdLog = zeros(size(sequenceBlockInitiationTimes));
    for behPrd = 1:size(behPeriods,1)
        seqBlockInitTimesBehPrdLog = seqBlockInitTimesBehPrdLog + (sequenceBlockInitiationTimes>behPeriods(behPrd,1) & sequenceBlockInitiationTimes<behPeriods(behPrd,2));
    end
    sequenceBlockInitiationTimes(logical(~seqBlockInitTimesBehPrdLog)) = [];
end

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

if term
    if ~multStart
        odorPresTermLog = cellfun(@(a)a<terminateTime, odorPresTime, 'uniformoutput', 0);
    elseif multStart
        tempOdorPresTermLog = cell(size(behPeriods,1),size(odorPresTime,2));
        for behPrd = 1:size(behPeriods,1)
            tempOdorPresTermLog(behPrd,:) = cellfun(@(a)a>behPeriods(behPrd,1) & a<behPeriods(behPrd,2), odorPresTime, 'uniformoutput', 0);
        end
        odorPresTermLog = cell(size(odorPresTime));
        for o = 1:length(odorPresTime)
            odorPresTermLog{o} = logical(sum(cell2mat(tempOdorPresTermLog(:,o)'),2));
        end
    end
    odorPresTime = cellfun(@(a,b)a(b), odorPresTime, odorPresTermLog, 'uniformoutput', 0);
end

if length(odorPresTime{1}) < length(odorPresTime{2})
%     plxSummary.Errors = {'More Bs than As'};
    error('More Bs than As')
%     plxData.Session = plxSession;
%     plxData.Summary = plxSummary;
    return 
end

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
        plxSummary.Errors = [plxSummary.Errors; {'More pokes initiated than ended'}];
        pokeInitiationTimes(end) = [];
    elseif isempty(pokeChanRep) && (length(pokeInitiationTimes) < length(pokeEndTimes))
        plxSummary.Errors = [plxSummary.Errors; {'More pokes ended than initiated'}];
        pokeEndTimes(1) = [];
    elseif ~isempty(pokeChanRep)
        if allPokes(pokeChanRep,2)==2
            pokeOutChanNumLog = pokeEndTimes==allPokes(pokeChanRep+1,1);
            plxSummary.Errors = [plxSummary.Errors; {'More pokes ended than initiated @' num2str(pokeEndTimes(pokeOutChanNumLog))}];
            pokeEndTimes(pokeOutChanNumLog) = [];
        elseif allPokes(pokeChanRep,2)==1
            pokeInitiationTimeLog = pokeInitiationTimes==allPokes(pokeChanRep,1);
            plxSummary.Errors = [plxSummary.Errors; {'More pokes initiated than ended @' num2str(pokeInitiationTimes(pokeInitiationTimeLog))}];
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

if term 
    if ~multStart
        pokesAfterTermLog = pokeInitiationTimes<terminateTime;
    elseif multStart
        tempPokesAfterTermLog = zeros(size(pokeInitiationTimes));
        for behPrd = 1:size(behPeriods,1)
            tempPokesAfterTermLog = tempPokesAfterTermLog + (pokeInitiationTimes>behPeriods(behPrd,1) & pokeInitiationTimes<behPeriods(behPrd,2));
        end
        pokesAfterTermLog = logical(tempPokesAfterTermLog);
    end
    pokeInitiationTimes = pokeInitiationTimes(pokesAfterTermLog);
    pokeDuration = pokeDuration(pokesAfterTermLog);
    interPokeInterval = interPokeInterval(pokesAfterTermLog);
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

if term
    if ~multStart
        frontRewardTermTimeLog = frontRewardTimes<terminateTime;
        backRewardTermTimeLog = backRewardTimes<terminateTime;
    elseif multStart
        tempFrontRewardTermTimeLog = zeros(size(frontRewardTimes));
        tempBackRewardTermTimeLog = zeros(size(backRewardTimes));
        for behPrd = 1:size(behPeriods,1)
            tempFrontRewardTermTimeLog = tempFrontRewardTermTimeLog + (frontRewardTimes>behPeriods(behPrd,1) & frontRewardTimes<behPeriods(behPrd,2));
            tempBackRewardTermTimeLog = tempBackRewardTermTimeLog + (backRewardTimes>behPeriods(behPrd,1) & backRewardTimes<behPeriods(behPrd,2));
        end
        frontRewardTermTimeLog = logical(tempFrontRewardTermTimeLog);
        backRewardTermTimeLog = logical(tempBackRewardTermTimeLog);
    end
    frontRewardTimes = frontRewardTimes(frontRewardTermTimeLog);
    backRewardTimes = backRewardTimes(backRewardTermTimeLog);
end

%% Position Data
% Convert the strobed channel into x-y coordinates using plexon
% functions
posChanNum = strcmp(channels, 'Strobed');
[~, ~, plxSummary.PositionVTmode, aniPosition] = plx_vt_interpret(plxStruct(posChanNum).ts, plxStruct(posChanNum).sv);

% clean position data of non-values
nonPositionLog = (aniPosition(:,2) + aniPosition(:,3))==0;
aniPosition(nonPositionLog,:) = [];

%% Identify and correct issues caused by asynchronous starts of Plexon and Matlab
if odorPresTime{1}(1)<pokeInitiationTimes(1) % This happens when there's a poke in the buffer prior to when Matlab starts
    % The solution is to start on the second sequence block
    newStartTime = sequenceBlockInitiationTimes(2);
    taintedOdorPresLog = cellfun(@(a)a<newStartTime, odorPresTime, 'uniformoutput', 0);
    odorPresTime = cellfun(@(a,b)a(~b), odorPresTime, taintedOdorPresLog, 'uniformoutput', 0);
    pokeInitiationTimes(pokeInitiationTimes<newStartTime) = [];
    sequenceBlockInitiationTimes(1) = [];
    plxSummary.Errors = [plxSummary.Errors; {'OdorPres1 < Poke1; Removed First Block'}];
end

firstSeqBlockStart = sequenceBlockInitiationTimes(1); % FSBS
if sum(pokeInitiationTimes<firstSeqBlockStart)>=1;
    preFSBSodorPres = cell2mat(cellfun(@(a)sum(a<firstSeqBlockStart), odorPresTime, 'uniformoutput', 0));
    if sum(preFSBSodorPres)==0 % This means there were no odors prior to first sequence block initiation, i.e. MATLAB was likely not running when the pokes were registered
        lastPokePreFSBS = pokeInitiationTimes(find(pokeInitiationTimes<firstSeqBlockStart, 1, 'last'));
        sequenceBlockInitiationTimes(1) = lastPokePreFSBS-0.1;
    elseif sum(preFSBSodorPres)>=1 % I believe this happens when Plexon was started after MATLAB
        sequenceBlockInitiationTimes = [pokeInitiationTimes(1)-0.1; sequenceBlockInitiationTimes];
    end
end

%% Clean up spurious odor times
odorTimesAll = cellfun(@(a,b)[a(:), ones(length(a),1).*b], odorPresTime(1:5), num2cell(1:5), 'uniformoutput',0);
odorTimesAllSorted = sortrows(cell2mat(odorTimesAll(:)));
odorTimesAllDiff = diff(odorTimesAllSorted(:,1));
spuriousOdors = sortrows([odorTimesAllSorted(([odorTimesAllDiff; 1]<1),:); odorTimesAllSorted(([1; odorTimesAllDiff]<1),:)]);
spuriousOdor1times = spuriousOdors((spuriousOdors(:,2)==1),1);

cleanOdorTimes = odorTimesAllSorted(~(ismember(odorTimesAllSorted(:,1),spuriousOdor1times)),:);
trialsWithSpuriousOdor1Activation = find(ismember(cleanOdorTimes(:,1),spuriousOdors(~(spuriousOdors(:,2)==1),1)));

sequenceLength = max(cleanOdorTimes(:,2));
plxSummary.SequenceLength = sequenceLength;

%% Fill in Trial Data
% Rather than going by blocks (as in the old code) this one just goes and
% looks around odor presentations and fills in everything around that.

ordPos = 1;
for trl = 1:length(cleanOdorTimes)
    %% Identify current trial's basic information
    curOdorTime = cleanOdorTimes(trl,1);
    if trl == length(cleanOdorTimes)
        nextOdorTime = curOdorTime + 2000;
        nextOdor = nan;
    else
        nextOdorTime = cleanOdorTimes(trl+1,1);
        nextOdor = cleanOdorTimes(trl+1,2);
    end
    plxSession(trl).OrdinalPosition = ordPos; % Ordinal position of the trial
    plxSession(trl).SequenceItem = cleanOdorTimes(trl,2); % Odor that was presented
    plxSession(trl).TranspositionDistance = ordPos - cleanOdorTimes(trl,2); % Relative distance of that item from InSeq
        
    % Determine which sequence block this trial belonged to
    sessionBlockStartTimePos = find(sequenceBlockInitiationTimes<curOdorTime,1,'last');
    if isempty(sessionBlockStartTimePos)
        error('Missing Sequence Block Start Time.... THIS SHOULD BE FIXED BY LINES 217-226... TELL GABE');
    end
    plxSession(trl).SessionBlockNumber = sessionBlockStartTimePos;
    
    %% Identify trial event times
    plxSession(trl).SessionBlockStartTime = sequenceBlockInitiationTimes(sessionBlockStartTimePos); % Time the trial block started
    plxSession(trl).ItemPresentationTime = curOdorTime; % The time the odor was presented
    frontReward = frontRewardTimes((frontRewardTimes>curOdorTime) & (frontRewardTimes<nextOdorTime));
    if isempty(frontReward)
        plxSession(trl).FrontRewardTime = nan;
    else
        plxSession(trl).FrontRewardTime = frontReward;
    end
    backReward = backRewardTimes((backRewardTimes>curOdorTime) & (backRewardTimes<nextOdorTime));
    if isempty(backReward)
        plxSession(trl).BackRewardTime = nan;
    else
        plxSession(trl).BackRewardTime = backReward;
    end
    
    %% Identify poke times
    preOdorPokesLog = pokeInitiationTimes <= curOdorTime;
    preOdorPokes = pokeInitiationTimes(preOdorPokesLog);
    plxSession(trl).OdorTrigPokeTime = max(preOdorPokes);
    odorTrigPokeDur = pokeDuration(pokeInitiationTimes==max(preOdorPokes));
    odorTrigPokeIPI = interPokeInterval(find(preOdorPokesLog,1,'last'));
    
    pokeOdorDelvOffset = plxSession(trl).ItemPresentationTime - plxSession(trl).OdorTrigPokeTime;
    
    if ~isempty(frontReward) && (plxSession(trl).TranspositionDistance == 0)
        trialEndTime = plxSession(trl).FrontRewardTime;
    elseif isempty(frontReward) && (plxSession(trl).TranspositionDistance == 0)
        trialEndTime = curOdorTime + odorTrigPokeDur;
    elseif isempty(frontReward) && ~(plxSession(trl).TranspositionDistance == 0)
        trialEndTime = curOdorTime + 1;
    end
    trialPokeLog = (pokeInitiationTimes>curOdorTime) & (pokeInitiationTimes<trialEndTime);
    trialPokes = pokeInitiationTimes(trialPokeLog);
    trialPokesDurations = pokeDuration(trialPokeLog);
    if ~(length(trialPokes)==length(trialPokesDurations))
        error('Go here this code should not. Tell gabe you should!');
    end
    trialInterPokeIntervals = interPokeInterval(trialPokeLog);
    if ~isempty(trialPokes)
        plxSession(trl).MultiOdorPokeLog = 1;
        plxSession(trl).OdorPokesDurations = [[max(preOdorPokes), odorTrigPokeDur, odorTrigPokeIPI]; [trialPokes, trialPokesDurations, trialInterPokeIntervals]];
        plxSession(trl).PokeDuration = sum(plxSession(trl).OdorPokesDurations(:,2));
        plxSession(trl).OdorPokeWithdrawTime = max(trialPokes+trialPokesDurations);
    else
        plxSession(trl).MultiOdorPokeLog = 0;
        plxSession(trl).OdorPokesDurations = [];
        plxSession(trl).PokeDuration = odorTrigPokeDur;
        plxSession(trl).OdorPokeWithdrawTime = max(preOdorPokes) + odorTrigPokeDur;
    end
    
    %% Determine Trial Outcomes
    plxSession(trl).Performance = ~isempty(frontReward); % The performance MATLAB determined
    switch ~isempty(frontReward)
        case 1
            if plxSession(trl).TranspositionDistance==0 && plxSession(trl).PokeDuration<1
                plxSummary.Errors = [plxSummary.Errors; {['Incorrect Duration Timing for trial ' num2str(trl) '. Was ' num2str(plxSession(trl).PokeDuration) ' on an InSeq Trial']}];
            elseif ~(plxSession(trl).TranspositionDistance==0) && plxSession(trl).PokeDuration>1
                plxSummary.Errors = [plxSummary.Errors; {['Incorrect Duration Timing for trial ' num2str(trl) '. Was ' num2str(plxSession(trl).PokeDuration) ' on an OutSeq Trial']}];
            end
            corRejNoFdBk = 0;
        case 0
            if ~(plxSession(trl).TranspositionDistance==0)
                if (plxSession(trl).PokeDuration < 1) && (nextOdor == ordPos + 1 || ordPos == sequenceLength)
                    corRejNoFdBk = 1;
                    plxSession(trl).Performance = 1;
                elseif (plxSession(trl).PokeDuration >1) && (nextOdor == ordPos +1)
                    corRejNoFdBk = 1;
                    plxSession(trl).Performance = 0;
                end
            else
                if plxSession(trl).TranspositionDistance==0 && plxSession(trl).PokeDuration>1
                    plxSummary.Errors = [plxSummary.Errors; {['Incorrect Duration Timing for trial ' num2str(trl) '. Was ' num2str(plxSession(trl).PokeDuration) ' on an InSeq Trial']}];
                elseif ~(plxSession(trl).TranspositionDistance==0) && plxSession(trl).PokeDuration<1
                    plxSummary.Errors = [plxSummary.Errors; {['Incorrect Duration Timing for trial ' num2str(trl) '. Was ' num2str(plxSession(trl).PokeDuration) ' on an OutSeq Trial']}];
                end
                corRejNoFdBk = 0;
            end
    end
                
    %% Identify trial logicals and increment/reset ordinal position tracker
    % Identify if trial had a spurious odorA activation
    if sum(trialsWithSpuriousOdor1Activation==trl)==1
        plxSession(trl).SpuriousOdorA = 1;
    else
        plxSession(trl).SpuriousOdorA = 0;
    end
    
    % Move along trial position
    if (isempty(frontReward) && ~corRejNoFdBk) || ordPos == sequenceLength
        ordPos = 1;
    else
        ordPos = ordPos + 1;
    end
end

%% Now go through and fill in the position data into trials
for trl = 1:length(plxSession)
    trlStart = plxSession(trl).OdorTrigPokeTime;
    if trl == length(plxSession)
        nextTrlStart = aniPosition(end,1);
    else
        nextTrlStart = plxSession(trl+1).OdorTrigPokeTime;
    end
    
    trialPosDataLog = aniPosition(:,1)>=trlStart & aniPosition(:,1)<nextTrlStart;
    curAniPosData = aniPosition(trialPosDataLog,:);
    plxSession(trl).PositionData = curAniPosData;
    aniInStartPos = (curAniPosData(:,2)>800 & curAniPosData(:,2)<900) & (curAniPosData(:,3)>310 & curAniPosData(:,3)<410);
    if sum(aniInStartPos)>=1
        plxSession(trl).AniReturnToOrigin = 1;
    else
        plxSession(trl).AniReturnToOrigin = 0;
    end
end

%% Compile Outputs
plxData.Raw = plxSession;
plxData.Summary = plxSummary;
