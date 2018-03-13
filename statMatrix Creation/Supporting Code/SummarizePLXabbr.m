function [plxData] = SummarizePLXabbr(sessionDataFile, sequenceLength, outInBuffer)
%% FindAndSummarizePLXfromMATfile
% Uses the string sessionDataFile that includes file name and directory.
% This function assumes files are organized in a consistent way, i.e. the
% .mat file being specified by sessionDataFile is located in a subfolder
% within an animal's data folder where there is a matching subfolder where
% the plx files are stored.
%
%%
global ssnFile % This is here so that if the script errors out I'll be able to tell which file it crapped out on
if nargin == 2 || nargin == 1
    outInBuffer = 0.1;
end
%% Identify plx path
if strcmp(sessionDataFile(end-2:end), 'plx')
    plxFile = sessionDataFile;
    fileLog = 1;
elseif strcmp(sessionDataFile(end-2:end), 'mat')
    ssnFilePathing = strsplit(sessionDataFile, '\');
    aniFile = ssnFilePathing{end};
    numBranches = length(ssnFilePathing);
    lengthBranches = cell2mat(cellfun(@(a)length(a), ssnFilePathing, 'uniformoutput', 0));
    aniDir = sessionDataFile(1:nansum(lengthBranches(1:numBranches-2))+(numBranches-2));
    aniPLXdir = [aniDir 'Plexon Files\'];
    aniPLXfile = [aniFile(1:end-16) '.plx'];
    
    plxDirFiles = dir(aniPLXdir);
    plxDirFileNames = {plxDirFiles.name};
    fileLog = cell2mat(cellfun(@(a)strcmp(a, aniPLXfile), plxDirFileNames , 'uniformoutput', 0));
    plxFile = [aniPLXdir aniPLXfile];
end
%% Create plxData Structure
% The structure plxData is a 1 x n structure where n = the number of trials
plxSession = struct('OrdinalPosition', [], 'SequenceItem', [], 'TranspositionDistance', [],...
    'SessionBlockNumber', [], 'Performance', [], 'SessionBlockStartTime', [],...
    'ItemPresentationTime', [], 'OdorTrigPokeTime', [], 'PokeDuration', [],...
    'OdorPokeWithdrawTime', [],...
    'OdorPokesDurations', [], 'RewardSignalTime', [], 'TrialEndTime', [],...
    'ErrorSignalTime', [], 'PositionData', [], 'AniReturnToOrigin', [],...
    'FrontRewardTime', [], 'BackRewardTime', [], 'MultiOdorPokeLog', []);
plxSummary.Errors = [];
%% Turn PLX file into structure
if nansum(fileLog) > 1
    error('Multiple PLX Files with same name, something wrong with file identification code and/or naming')
elseif nansum(fileLog) == 1
    % Compile Plexon Channels into a single structure
    % This may be a but unnecessary given the hard coding I do immediately
    % after this but a) it's easier to work with this structure and b)
    % automating the channel data this way allows for a subfunction to work
    % with new channel definitions should the event channels change...
    % which they probably will when the rig is upgraded to have 10 odors
    [numChans, chanNames] = plx_event_names(plxFile);
    plxStruct = struct('Channel', [], 'n', [], 'ts', [], 'sv', []);
    plxSsn = [];
    for chan = 1:numChans
        curChan = chanNames(chan,:);
        intVals = double(curChan);
        valLim = find(~(intVals==0), 1, 'last');
        plxStruct(chan).Channel = curChan(1:valLim);
        [plxStruct(chan).n, plxStruct(chan).ts, plxStruct(chan).sv] = plx_event_ts(plxFile, curChan);
        plxSsn = [plxSsn; [plxStruct(chan).ts, ones(size(plxStruct(chan).ts)).*chan]];
    end
    plx_close(plxFile);
    channels = {plxStruct.Channel};
    strobeVal = {plxStruct.sv};
    
    %% Check for Activity on the Terminate Channel
    % Read the terminate channel and determine if a terminate command was
    % issued anytime during the session. If so then truncate all the other
    % behavioral variables to only the time points that occurred prior to
    % when the command was issued
    terminateChannelLog = strcmp(channels, 'Terminate');
    if ~(sum(terminateChannelLog)==0) && ~(plxStruct(terminateChannelLog).n == 0)
        term = 1;
        plxSummary.Terminate = 1;
        terminateTime = plxStruct(terminateChannelLog).ts(find(plxStruct(terminateChannelLog).ts>0, 1, 'first'));
    else
        term = 0;
        plxSummary.Terminate = 0;
    end
    
    %% Pull Information from Buzzer channel (Sequence Block, Trial and Error Identifiers)
    % The easiest to parse the session is using the buzzer channel as it
    % has unique ways/times of activating during the trial.
    % There are two modes of activation:
    % 1) Double Buzzer: Occurs when the animal is in the back of the maze
    %       to start a new sequence
    % 2) Single Buzzer: Occurs in between trials as an indicator that the
    %       animal can poke to start a new trial (short activation of the
    %       buzzer) as well as that the animal made an error in the
    %       sequence (long activation of the buzzer).
    % It is impossible to distinguish the two kinds of single buzzer
    % activation apart from the context in which they occur. Specifically
    % it's the activity on the reward channels that disambiguates the short
    % vs. long activation:
    %
    % Short Activation (trial): Only occurs following a correct and is only
    %   followed by a double when there are front and back water rewards in
    %   between the two. The reward signal is an unreliable indicator of
    %   reward presentation as it doesn't activate on correct rejection
    %   (correct OutSeq) trials.
    % Long Activation (error): Follows an incorrect trial and is followed
    %   by a double buzzer or the end of the session. There will never be
    %   signal on either reward channel between it and the double buzzer.
    %
    % The double buzzer is easily distinguished from either short
    % activation case because it consists of two activations on the
    % buzzer channel within 0.1s of eachother. 0.12s is used here just to
    % be on the safe side.
    
    buzzerChanNum = strcmp(channels, 'Buzzer');
    interSolenoidInterval1 = diff([plxStruct(buzzerChanNum).ts; 0]); % Inter-Solenoid Interval to identify first double buzz
    interSolenoidInterval2 = diff([0; plxStruct(buzzerChanNum).ts]); % Inter-Solenoid Interval to identify second double buzz
    sequenceBlockInitiationSpots = interSolenoidInterval1<=0.2 & interSolenoidInterval1>0; % Identify the first of the two double buzzes
    sequenceBlockInitiationTimes = plxStruct(buzzerChanNum).ts(sequenceBlockInitiationSpots);
    nonDoubleBuzzBuzzer = plxStruct(buzzerChanNum).ts(interSolenoidInterval1>0.2 & interSolenoidInterval2>0.2);

    clear interSolenoidInterval sequenceBlockInitiationSpots
    
    if term
        sequenceBlockInitiationTimes(sequenceBlockInitiationTimes>terminateTime) = [];
        nonDoubleBuzzBuzzer(nonDoubleBuzzBuzzer>terminateTime) = [];
    end
    %% Identify Beep (Reward Signal) Times
    beepChanNum = strcmp(channels, 'Tone');
    beepTimes = plxStruct(beepChanNum).ts;
    
    if term 
        beepTimes(beepTimes>terminateTime) = [];
    end
    %% Identify Trial Odors
    % Pull out odor presentation times
    odorPresTime = cell(1,10); % ASSUMES 10 ODORS; CHANGE IF THIS ASSUMPTION IS NO LONGER CORRECT
    
    odorAchanNum = strcmp(channels, 'Odor A');
    odorAchanTS = plxStruct(odorAchanNum).ts;
    odorAlog = plxStruct(odorAchanNum).sv==0;
    odorVlog = plxStruct(odorAchanNum).sv==1;
    odorPresTime{1} = odorAchanTS(odorAlog);
    odorPresTime{6} = odorAchanTS(odorVlog);
    
    odorBchanNum = strcmp(channels, 'Odor B');
    odorBchanTS = plxStruct(odorBchanNum).ts;
    odorBlog = plxStruct(odorBchanNum).sv==0;
    odorWlog = plxStruct(odorBchanNum).sv==1;
    odorPresTime{2} = odorBchanTS(odorBlog);
    odorPresTime{7} = odorBchanTS(odorWlog);
    
    odorCchanNum = strcmp(channels, 'Odor C');
    if ~sum(odorCchanNum)==0
        odorCchanTS = plxStruct(odorCchanNum).ts;
        odorClog = plxStruct(odorCchanNum).sv==0;
        odorXlog = plxStruct(odorCchanNum).sv==1;
        odorPresTime{3} = odorCchanTS(odorClog);
        odorPresTime{8} = odorCchanTS(odorXlog);
    end
    
    odorDchanNum = strcmp(channels, 'Odor D');
    if ~sum(odorDchanNum)==0
        odorDchanTS = plxStruct(odorDchanNum).ts;
        odorDlog = plxStruct(odorDchanNum).sv==0;
        odorYlog = plxStruct(odorDchanNum).sv==1;
        odorPresTime{4} = odorDchanTS(odorDlog);
        odorPresTime{9} = odorDchanTS(odorYlog);
    end
    
    odorEchanNum = strcmp(channels, 'Odor E');
    if ~sum(odorEchanNum)==0
        odorEchanTS = plxStruct(odorEchanNum).ts;
        odorElog = plxStruct(odorEchanNum).sv==0;
        odorZlog = plxStruct(odorEchanNum).sv==1;
        odorPresTime{5} = odorEchanTS(odorElog);
        odorPresTime{10} = odorEchanTS(odorZlog);
    end
    
    if term
        odorPresTermLog = cellfun(@(a)a>terminateTime, odorPresTime, 'uniformoutput', 0);
        odorPresTime = cellfun(@(a,b)a(~b), odorPresTime, odorPresTermLog, 'uniformoutput', 0);
    end
    
    if length(odorPresTime{1}) < length(odorPresTime{2})
        plxSummary.Errors = {'More Bs than As'};
%         error('More Bs than As')
        plxData.Session = plxSession;
        plxData.Summary = plxSummary;
        return
    end
    
    if ~exist('sequenceLength', 'var')
        sequenceLength = sum(cellfun(@(a)~isempty(a), odorPresTime));
    end
    plxSummary.SequenceLength = sequenceLength;

    %% Identify Poke Times
    % The key information used here is when the poke was initiated and how
    % long it was held. When the poke was broken is less important and
    % easily derived from (initiation time) + (hold time).
    pokeInChanNum = strcmp(channels, 'Photocell IN');
    pokeInitiationTimes = plxStruct(pokeInChanNum).ts;
    pokeOutChanNum = strcmp(channels, 'Photocell OUT') | strcmp(channels, 'Photocell Out');
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
    pokePairs = [pokeInitiationTimes pokeEndTimes];
    
    interPokeInterval = diff([[0; pokeEndTimes], [pokeInitiationTimes; 0]],1,2);
    interPokeInterval(1) = []; % Remove first inter-poke interval because it should refer to period between first and second poke, not period from start to end of first poke
    interPokeInterval(end) = nan; % NaN out the final entry because there is no poke in after the final poke out
    if nansum(interPokeInterval<0)>0
        error('Inter Poke Intervals are not all positive, something is wrong');
    end
    
    if term
        pokeInitiationTimes(pokeInitiationTimes>terminateTime) = [];
        pokeDuration(pokeInitiationTimes>terminateTime) = [];
        interPokeInterval(pokeInitiationTimes>terminateTime) = [];
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
    waterAchanNum = strcmp(channels, 'Water A (front)');
    interSolenoidInterval = diff([0; plxStruct(waterAchanNum).ts]); % Identify how long it was in between solenoid activations. The zero here is just so that in the following steps 1 is the first spot.
    firstSolenoidActivationTime = interSolenoidInterval>=1; % The water solenoid is activated 4x times when reward is triggered with <0.5s in between activations. This line finds the times the solenoid when >=1s between activations, i.e. in between reward presentations
    frontRewardTimes = plxStruct(waterAchanNum).ts(firstSolenoidActivationTime);
    clear interSolenoidInterval firstSolenoidActivationTime
    
    % Rear Reward
    waterBchanNum = strcmp(channels, 'Water B (back)');
    interSolenoidInterval = diff([0; plxStruct(waterBchanNum).ts]); % Identify how long it was in between solenoid activations. The zero here is just so that in the following steps 1 is the first spot.
    firstSolenoidActivationTime = interSolenoidInterval>=1; % The water solenoid is activated 4x times when reward is triggered with <0.5s in between activations. This line finds the times the solenoid when >=1s between activations, i.e. in between reward presentations
    backRewardTimes = plxStruct(waterBchanNum).ts(firstSolenoidActivationTime);
    clear interSolenoidInterval firstSolenoidActivationTime
    
    if term
        frontRewardTimes(frontRewardTimes>terminateTime) = [];
        backRewardTimes(frontRewardTimes>terminateTime) = [];
    end
    
    %% Identify Trial Performance
    corrTrialChanNum = strcmp(channels, 'Correct Poke') | strcmp(channels, 'Keyboard7');
    correctTrialsTimes = plxStruct(corrTrialChanNum).ts;
    
    inCorrTrialChanNum = strcmp(channels, 'Incorrect Poke') | strcmp(channels, 'Keyboard8');
    incorrectTrialsTimes = plxStruct(inCorrTrialChanNum).ts;
    
    if term
        correctTrialsTimes(correctTrialsTimes>terminateTime) = [];
        incorrectTrialsTimes(incorrectTrialsTimes>terminateTime) = [];
    end
    
    trialOutcomes = sort([correctTrialsTimes;incorrectTrialsTimes]);
    if sum(pokeInitiationTimes(1) > trialOutcomes) >= 1 % If a trial outcome is registered before a poke is registered that means there was a poke already in the buffer and matlab saw that poke in but not the accompanying poke out. This happens when the nose port is cleaned between sessions and the buffer is not cleared before running the next session.
        miscreant = trialOutcomes(pokeInitiationTimes(1) > trialOutcomes);
        correctTrialsTimes(correctTrialsTimes==miscreant) = [];
        incorrectTrialsTimes(incorrectTrialsTimes==miscreant) = [];
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
    if sum(pokeInitiationTimes<firstSeqBlockStart)>=1
        preFSBSodorPres = cell2mat(cellfun(@(a)sum(a<firstSeqBlockStart), odorPresTime, 'uniformoutput', 0));
        if sum(preFSBSodorPres)==0 % This means there were no odors prior to first sequence block initiation, i.e. MATLAB was likely not running when the pokes were registered
            lastPokePreFSBS = pokeInitiationTimes(find(pokeInitiationTimes<firstSeqBlockStart, 1, 'last'));
            sequenceBlockInitiationTimes(1) = lastPokePreFSBS-0.1;
        elseif sum(preFSBSodorPres)>=1 % I believe this happens when Plexon was started after MATLAB
            sequenceBlockInitiationTimes = [pokeInitiationTimes(1)-0.1; sequenceBlockInitiationTimes];
        end
    end               
    
    
    
    %% Isolate Sequence Blocks & Identify Trials
    % The buzzer channel is the best way to parse sequence blocks.
    % Unfortunately it is an unreliable marker of trial start times. In
    % theory it signals start of the next trial however there are instances
    % where the time codes are messed up.
    % 
    % Specifically, if the ITI elapses and the animal pokes after the ITI
    % but prior to the buzzer that poke is registered in plexon prior to
    % buzzer but only read into MATLAB after the buzzer happens because of
    % the buffer. MATLAB then reads this poke and triggers odor
    % presentation quickly after the buzzer. 
    % 
    % Therefore, the most reliable trial marker is the odor presentation.
    %
    % Initially I had hoped to break the trial into pre-odor, odor and
    % post-odor periods. Due to the quirks of working with the plexon
    % events buffer (and how the code deals with that) it's best to just
    % focus on what we need. 1) When do they poke in, 2) How long is that
    % poke, 3) Do they poke out? 4) If they poke out, do they poke back in
    % within the poke buffer (100ms) period? 5) If so, how long is the poke
    % including the poke buffer?

    trl = 1;
    sb = 1;
    sbit = length(sequenceBlockInitiationTimes);
    while sb <= sbit
        emptyBlock = 0;
        %% Identify Sequence Block Boundaries
        seqBlockStart = sequenceBlockInitiationTimes(sb);
        
        if sb == length(sequenceBlockInitiationTimes)
            seqBlockEnd = sequenceBlockInitiationTimes(sb) + 2400; % If you're working on the final sequence block just assume the block ends 40min after when it starts. If the animal was just sitting there for 40min then the session probably should have been terminated anyway.
        else
            seqBlockEnd = sequenceBlockInitiationTimes(sb+1)-0.001;
        end
        %% Identify Odors and Odor Presentation Times during the block
        % Note, items remain in cells because multiple instances of the
        % same odor can be presented during the same block (skip/repeat).
        % If odors are not presented during the block their column will
        % contain an empty cell.
        seqBlockOdorPresTimes = cellfun(@(a)a(a>seqBlockStart & a<seqBlockEnd), odorPresTime, 'uniformoutput', 0);
        seqBlockMptyLog = cell2mat(cellfun(@(a)isempty(a), seqBlockOdorPresTimes, 'uniformoutput',0));
        seqBlockOdorPresTimes(seqBlockMptyLog) = {[]};
        for si = 1:length(seqBlockOdorPresTimes)
            if ~isempty(seqBlockOdorPresTimes{si})
                seqBlockOdorPresTimes{si} = [seqBlockOdorPresTimes{si} ones(size(seqBlockOdorPresTimes{si})).*si];
            end
        end
        seqBlockOdors = sortrows(cell2mat(seqBlockOdorPresTimes'));
        if size(seqBlockOdors,1)>sequenceLength
            if nansum(diff(sort(seqBlockOdors(:,1)))<0.1) == sequenceLength % This happens when the plexon file isn't stopped before the odor system was depressurized and the command was captured by the recording
                plxSummary.Errors = [plxSummary.Errors; {'Captured Pressure Release'}];
                for o = 1:length(sequenceLength)
                    seqBlockOdorPresTimes{o}(end,:) = [];
                    if isempty(seqBlockOdorPresTimes{o})
                        seqBlockOdorPresTimes{o} = [];
                    end
                end
                seqBlockOdors = sortrows(cell2mat(seqBlockOdorPresTimes'));
            else
                sequenceBlockInitiationTimes = sort([sequenceBlockInitiationTimes; seqBlockOdors(sequenceLength+1,1)-0.0001]);
                seqBlockOdors = seqBlockOdors(1:sequenceLength,:);
            end
        end
        if isempty(seqBlockOdors)
            if emptyBlock == 0
                emptyBlockTime = sequenceBlockInitiationTimes(sb);
                sequenceBlockInitiationTimes(sb) = [];
                emptyBlock = 1;
                sbit = sbit - 1;
            end
            plxSummary.Errors = [plxSummary.Errors ; {['Empty Block @' num2str(emptyBlockTime) '; No Odors']}];
%             if sb == length(sequenceBlockInitiationTimes)
%                 break
% %                 error('Why you get here code? Go Home code, no one wants you here!')
%             else
%                 %             break
%                 error('No Odors!!! THIS SHOULD BE FIXED!!!!')
%             end
        end
        %% Identify Single Buzzer Activation (Trial & Error Signals) Times during the block
        seqBlockBuzz = nonDoubleBuzzBuzzer(nonDoubleBuzzBuzzer>seqBlockStart & nonDoubleBuzzBuzzer<seqBlockEnd);
        
        %% Identify Beep Activation (Reward Signal) and Reward Presentation Times during the block
        seqBlockRewardSignal = beepTimes(beepTimes>seqBlockStart & beepTimes<seqBlockEnd);
        seqBlockFrontRewardTimes = frontRewardTimes(frontRewardTimes>seqBlockStart & frontRewardTimes<seqBlockEnd);
        seqBlockBackRewardTimes = backRewardTimes(backRewardTimes>seqBlockStart & backRewardTimes<seqBlockEnd);
        
        %% Identify Poke Times and Durations during the block
        % Note that the hold duration selection is based on the logic of
        % when the poke was initiated.
        seqBlockPokeTimes = pokeInitiationTimes(pokeInitiationTimes>=seqBlockStart & pokeInitiationTimes<seqBlockEnd);
        seqBlockPokeDurations = pokeDuration(pokeInitiationTimes>seqBlockStart & pokeInitiationTimes<seqBlockEnd);
        seqBlockInterPokeIntervals = interPokeInterval(pokeInitiationTimes>seqBlockStart & pokeInitiationTimes<seqBlockEnd);
        if isempty(seqBlockPokeTimes)
            if emptyBlock == 0
                emptyBlockTime = sequenceBlockInitiationTimes(sb);
                sequenceBlockInitiationTimes(sb) = [];
                emptyBlock = 1;
                sbit = sbit - 1;
            end
            plxSummary.Errors = [plxSummary.Errors; {['Empty Block @' num2str(emptyBlockTime) '; No Pokes']}];
        end
        %% Identify Trial Results
        % What was the trial identified as?
        seqBlockCorrectTrialsTimes = correctTrialsTimes(correctTrialsTimes>seqBlockStart & correctTrialsTimes<seqBlockEnd);
        if isempty(seqBlockCorrectTrialsTimes)
            seqBlockCorrectTrialsTimes = [];
        end
        seqBlockIncorrectTrialsTimes = incorrectTrialsTimes(incorrectTrialsTimes>seqBlockStart & incorrectTrialsTimes<seqBlockEnd);
        if isempty(seqBlockIncorrectTrialsTimes)
            seqBlockIncorrectTrialsTimes = [];
        end
        seqBlockPerf = sortrows([[seqBlockCorrectTrialsTimes, ones(size(seqBlockCorrectTrialsTimes))];[seqBlockIncorrectTrialsTimes, zeros(size(seqBlockIncorrectTrialsTimes))]]);
        
        if size(seqBlockPerf,1) < size(seqBlockOdors,1)
            % Find the missing trial performance value
            tempTrlcounter = 1;
            while tempTrlcounter <= size(seqBlockOdors,1)                
                curPerfTime = seqBlockPerf(tempTrlcounter,1);
                if sum(curPerfTime < seqBlockOdors(:,1)) == size(seqBlockOdors,1)-tempTrlcounter
                elseif sum(curPerfTime < seqBlockOdors(:,1)) > size(seqBlockOdors,1)-tempTrlcounter
                    error('WHAT FUCKERY IS THIS');
                elseif sum(curPerfTime < seqBlockOdors(:,1)) < size(seqBlockOdors,1)-tempTrlcounter
                    if size(seqBlockFrontRewardTimes,1) == size(seqBlockOdors,1)
                        tempSeqBlockLaterPerf = seqBlockPerf(tempTrlcounter:end,:);
                        seqBlockPerf = seqBlockPerf(1:tempTrlcounter-1,:);
                        seqBlockPerf(tempTrlcounter,1) = seqBlockFrontRewardTimes(tempTrlcounter,1);
                        seqBlockPerf(tempTrlcounter,2) = 1;
                        seqBlockPerf = [seqBlockPerf; tempSeqBlockLaterPerf]; %#ok<AGROW>
                        break
                    elseif (size(seqBlockFrontRewardTimes,1) == size(seqBlockPerf,1))  && (seqBlockPerf(end,2) == 0)
                        tempSeqBlockLaterPerf = seqBlockPerf(tempTrlcounter:end,:);
                        seqBlockPerf = seqBlockPerf(1:tempTrlcounter-1,:);
                        seqBlockPerf(tempTrlcounter,1) = seqBlockFrontRewardTimes(tempTrlcounter,1);
                        seqBlockPerf(tempTrlcounter,2) = 1;
                        seqBlockPerf = [seqBlockPerf; tempSeqBlockLaterPerf]; %#ok<AGROW>
                    end
                else
                    error('Fuck this code');
                end
                if tempTrlcounter == size(seqBlockPerf,1)
                    if size(seqBlockFrontRewardTimes,1) == size(seqBlockOdors,1)
                        seqBlockPerf(tempTrlcounter+1,1) = seqBlockFrontRewardTimes(tempTrlcounter+1,1);
                        seqBlockPerf(tempTrlcounter+1,2) = 1;
                        break
                    elseif size(seqBlockPerf,1) == size(seqBlockOdors,1)
                        break
                    else
                        error('I hate this code');
                    end
                else
                    tempTrlcounter = tempTrlcounter + 1;
                end
            end
        end
                
        
        if isempty(seqBlockPerf)
            if emptyBlock == 0
                emptyBlockTime = sequenceBlockInitiationTimes(sb);
                sequenceBlockInitiationTimes(sb) = [];
                emptyBlock = 1;
                sbit = sbit - 1;
            end
            plxSummary.Errors = [plxSummary.Errors; {['Empty Block @' num2str(emptyBlockTime) '; No Performance']}];
%             if sb == length(sequenceBlockInitiationTimes)
%                 break
%             else
%                 error('No Performance!')
%             end
        end
        %% Break the Block into Trials and Analyze them Individually
        % For each odor presentation do the following steps to get poke
        % data:
        % 1) Identify odor presentation time
        % 2) Identify the poke that preceded odor presentation
        % 3) Identify when the trial ended (reward or error signal
        % presentation time)
        % 4) Identify whether additional pokes occurred during the odor
        % presentation that may have inter-poke intervals within the buffer
        % period.
        % 5) Determine poke duration and trial identification
        %         if sb==2; break; end
        if ~emptyBlock
            for t = 1:size(seqBlockOdors,1)
                %% Fill in Trial Identifiers
                plxSession(trl).SessionBlockNumber = sb;
                plxSession(trl).SessionBlockStartTime = seqBlockStart;
                plxSession(trl).OrdinalPosition = t;
                plxSession(trl).SequenceItem = seqBlockOdors(t,2);
                plxSession(trl).TranspositionDistance = t - plxSession(trl).SequenceItem;
                
                %% Identify Odor Presentation Time
                plxSession(trl).ItemPresentationTime = seqBlockOdors(t,1);
                
                %% Identify Poke Nearest to Odor Presentation
                preOdorPokes = seqBlockPokeTimes(seqBlockPokeTimes <= plxSession(trl).ItemPresentationTime);
                if isempty(preOdorPokes)
                    if ~isempty(pokeInitiationTimes(pokeInitiationTimes<=plxSession(trl).ItemPresentationTime))
                        lastPokeBeforeDoubleBuzzLog = pokeInitiationTimes==max(pokeInitiationTimes(pokeInitiationTimes<=plxSession(trl).ItemPresentationTime & pokeInitiationTimes>=plxSession(trl-1).ItemPresentationTime));
                        seqBlockPokeTimes = [pokeInitiationTimes(lastPokeBeforeDoubleBuzzLog); seqBlockPokeTimes];
                        seqBlockPokeDurations = [pokeDuration(lastPokeBeforeDoubleBuzzLog); seqBlockPokeDurations];
                        seqBlockInterPokeIntervals = [interPokeInterval(lastPokeBeforeDoubleBuzzLog); seqBlockInterPokeIntervals];
                        if t==1
                            preOdorPokes = seqBlockPokeTimes(seqBlockPokeTimes <= plxSession(trl).ItemPresentationTime);
                            if max(preOdorPokes) == plxSession(trl-1).ItemPresentationTime
                                error('Misc Err A.... FIX');
                            end
                        else
                            error('Misc Err B.... FIX');
                        end
                    elseif trl==1 && isempty(pokeInitiationTimes(pokeInitiationTimes<=plxSession(trl).ItemPresentationTime))
                        seqBlockPokeTimes = [(plxSession(trl).ItemPresentationTime-0.05); seqBlockPokeTimes];
                        seqBlockPokeDurations = [1; seqBlockPokeDurations];
                        seqBlockInterPokeIntervals = [1; seqBlockInterPokeIntervals];
                        preOdorPokes = seqBlockPokeTimes(seqBlockPokeTimes <= plxSession(trl).ItemPresentationTime);
                    end
                    if isempty(preOdorPokes)
                        error('Immaculate odor!')
                    end
                end
                plxSession(trl).OdorTrigPokeTime = max(preOdorPokes);
                
                %% Identify when odor trial ended
                if t == size(seqBlockOdors,1)
                    if sb == length(sequenceBlockInitiationTimes)
                        nextTrialStart = sequenceBlockInitiationTimes(sb)+3000;
                    else
                        nextTrialStart = sequenceBlockInitiationTimes(sb+1)-0.001;
                    end
                else
                    nextTrialStart = seqBlockOdors(t+1,1);
                end
                trialPerfLog = seqBlockPerf(:,1)>seqBlockOdors(t,1) & seqBlockPerf(:,1)<nextTrialStart;
                if nansum(trialPerfLog)>=2
                    error('Multiple Performance Indicators!')
                end
                plxSession(trl).Performance = seqBlockPerf(trialPerfLog,2);
                
                if plxSession(trl).Performance
                    trialEndTime = seqBlockFrontRewardTimes(seqBlockFrontRewardTimes>seqBlockOdors(t,1) & seqBlockFrontRewardTimes<nextTrialStart);
                    plxSession(trl).FrontRewardTime = trialEndTime;
                    backRewardTime = seqBlockBackRewardTimes(seqBlockBackRewardTimes>seqBlockOdors(t,1) & seqBlockBackRewardTimes<nextTrialStart);
                    if isempty(backRewardTime)
                        plxSession(trl).BackRewardTime = nan;
                    else
                        plxSession(trl).BackRewardTime = backRewardTime;
                    end
                    rewardSignalTime = seqBlockRewardSignal(seqBlockRewardSignal>seqBlockOdors(t,1) & seqBlockRewardSignal<nextTrialStart);
                    if isempty(rewardSignalTime)
                        plxSession(trl).RewardSignalTime = nan;
                    else
                        plxSession(trl).RewardSignalTime = rewardSignalTime;
                    end
                    plxSession(trl).ErrorSignalTime = nan;
                    plxSession(trl).TrialEndTime = trialEndTime;
                else
                    trialEndTime = seqBlockBuzz(seqBlockBuzz>seqBlockOdors(t,1) & seqBlockBuzz<nextTrialStart);
                    if length(trialEndTime)>1
                        trialEndTime = min(trialEndTime);
                        plxSession(trl).ErrorSignalTime = trialEndTime;
                    elseif isempty(trialEndTime)
                        if sb==length(sequenceBlockInitiationTimes)
                            trialEndTime = seqBlockEnd;
                            plxSession(trl).ErrorSignalTime = nan;
                        else
                            error('No-Buzz! No-Buzz!');
                        end
                    else
                        plxSession(trl).ErrorSignalTime = trialEndTime;
                    end
                    plxSession(trl).FrontRewardTime = nan;
                    plxSession(trl).BackRewardTime = nan;
                    plxSession(trl).RewardSignalTime = nan;
                    plxSession(trl).TrialEndTime = trialEndTime;
                end
                
                %% Identify Trial Pokes
                trialPokes = seqBlockPokeTimes(seqBlockPokeTimes>=plxSession(trl).OdorTrigPokeTime & seqBlockPokeTimes<trialEndTime);
                if isempty(trialPokes)
                    error('No Trial Pokes!')
                end
                trialPokesDurations = seqBlockPokeDurations(seqBlockPokeTimes>=plxSession(trl).OdorTrigPokeTime & seqBlockPokeTimes<trialEndTime);
                trialInterPokeIntervals = seqBlockInterPokeIntervals(seqBlockPokeTimes>=plxSession(trl).OdorTrigPokeTime & seqBlockPokeTimes<trialEndTime);
                if nansum(trialInterPokeIntervals>0.1)==0
                    plxSession(trl).PokeDuration = nansum(trialPokesDurations);
                    if length(trialPokes)>1
                        plxSession(trl).MultiOdorPokeLog = 1;
                    else
                        plxSession(trl).MultiOdorPokeLog = 0;
                    end
                else
                    longIPIspot = find(trialInterPokeIntervals>outInBuffer);
                    plxSession(trl).PokeDuration = nansum(trialPokesDurations(1:longIPIspot));
                    if length(trialPokes)>1
                        plxSession(trl).MultiOdorPokeLog = 1;
                    else
                        plxSession(trl).MultiOdorPokeLog = 0;
                    end
                end
                plxSession(trl).OdorPokeWithdrawTime = pokePairs(pokePairs(:,1)==trialPokes(end),2);
                plxSession(trl).OdorPokesDurations = trialPokesDurations;
                %% Fill in Performance Matrices
                if plxSession(trl).Performance
                    if plxSession(trl).TranspositionDistance==0 % True Correct
                        if plxSession(trl).PokeDuration < 1
                            plxSummary.Errors = [plxSummary.Errors; {['CorInSeq Hold<1 on Trial ' num2str(trl)]}];
                        else
                        end
                    elseif ~plxSession(trl).TranspositionDistance==0 % Correct Rejection
                        if  plxSession(trl).PokeDuration >= 1
                            plxSummary.Errors = [plxSummary.Errors; {['CorOutSeq Hold>1 on Trial ' num2str(trl)]}];
                        else
                        end
                    end
                elseif ~plxSession(trl).Performance
                    if plxSession(trl).TranspositionDistance==0 % False Negative
                        if plxSession(trl).PokeDuration >= 1
                            plxSummary.Errors = [plxSummary.Errors; {['InCorInSeq Hold>1 on Trial ' num2str(trl)]}];
                        else
                        end
                    elseif ~plxSession(trl).TranspositionDistance==0 % False Positive
                        if plxSession(trl).PokeDuration < 1
                            plxSummary.Errors = [plxSummary.Errors; {['InCorOutSeq Hold<1 on Trial ' num2str(trl)]}];
                        else
                        end
                    end
                end                
                %% Move to the next trial
                trl = trl+1;
            end
            sb = sb+1;
        else
        end
    end
    if isnan(plxSession(end).ErrorSignalTime) && isnan(plxSession(end).RewardSignalTime)
        plxSession(end) = [];
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
   

    
else
%     error(['No PLX file found for ' ssnFile '. Check your rig configuration and file name and or directory organization']);
    display(['No PLX file found for ' ssnFile '. Check your rig configuration and file name and or directory organization']);
    plxSummary.Errors = {'No PLX file'};
end
%% Compile Outputs
plxData.Raw = plxSession;
plxData.Summary = plxSummary;