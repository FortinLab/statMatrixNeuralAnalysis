function [plxData] = SummarizePLXabbr_BOS3(plxFile)
%% SummarizePLXabbr_BOS3
% Uses the string sessionDataFile that includes file name and directory.
% This is a version of SummarizePLXabbr.m designed to work with the data
% collected in the Eichenbaum lab that uses different event channels
%
% Version 2 created after it was recognized the poke durations were not
% accurate in version 1.
% Version 3 created after it was recognized that version 2 was not
% adequately capturing instances of buffer errors at poke withdrawal

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

%% Run through the session
for trl = 1:length(plxSession)
    clear trlPokesMtrx
    % Extract relevant trial variables
    curTargDur = txtTargDur(trl);
    curTxtBuffer = txtBuffer(trl)/1000;
    
    % Identify the start of the trial
    plxSession(trl).ItemPresentationTime = odorTimesAllSorted(trl,1);
    plxSession(trl).OdorTrigPokeTime =  pokeInitiationTimes(find(pokeInitiationTimes<odorTimesAllSorted(trl,1), 1, 'last'));
    plxSession(trl).OdorDeliveryDelay = plxSession(trl).ItemPresentationTime - plxSession(trl).OdorTrigPokeTime;
    
    % Now pull out reward instances & identify potential trial pokes (all
    % poke initiations between trial N and trial N+1)
    if trl == length(plxSession)
        frontRwdTime = frontRewardTimes(find(frontRewardTimes>odorTimesAllSorted(trl,1),1,'first'));
        rearRwdTime = backRewardTimes(find(backRewardTimes>odorTimesAllSorted(trl,1),1,'first'));
        allTrlPokesLog = pokeInitiationTimes>=plxSession(trl).OdorTrigPokeTime;
    else
        frontRwdTime = frontRewardTimes(find((frontRewardTimes>odorTimesAllSorted(trl,1)) & (frontRewardTimes<odorTimesAllSorted(trl+1,1)),1,'first'));
        rearRwdTime = backRewardTimes(find((backRewardTimes>odorTimesAllSorted(trl,1)) & (backRewardTimes<odorTimesAllSorted(trl+1,1)),1,'first'));
        allTrlPokesLog = (pokeInitiationTimes>=plxSession(trl).OdorTrigPokeTime) &...
            (pokeInitiationTimes<odorTimesAllSorted(trl+1,1));
    end
    if ~isempty(frontRwdTime)
        plxSession(trl).FrontRewardTime = frontRwdTime;
    end
    if ~isempty(rearRwdTime)
        plxSession(trl).BackRewardTime = rearRwdTime;
    end
        
    % Create the trial poke matrix, depicting the observed poking behavior
    % on the trial.
    % Matrix is organized as:
    % Column 1: Poke In Time
    % Column 2: Poke Out Times
    % Column 3: Poke In Times Relative to Trial Initiating Poke
    % Column 4: Poke Out Times Relative to Trial Initiating Poke
    % Column 5: Inter-Poke Interval (time till next poke in)
    trlPokesMtrx = [pokeInitiationTimes(allTrlPokesLog),...
        pokeEndTimes(allTrlPokesLog),...
        pokeInitiationTimes(allTrlPokesLog)-plxSession(trl).OdorTrigPokeTime,...
        pokeEndTimes(allTrlPokesLog)-plxSession(trl).OdorTrigPokeTime,...
        interPokeInterval(allTrlPokesLog)];
    % The way these timestamps are pulled, the last timestamp is always the
    % poke initiating the subsequent trial (unless it's the final trial in
    % the session) so remove it from this trial's poke log.
    if trl ~= length(plxSession)
        trlPokesMtrx(end,:) = [];
    end
    plxSession(trl).TrialPokesDurations = trlPokesMtrx;
    
    %% For Trials with only a single poke
    if size(trlPokesMtrx,1)==1
        plxSession(trl).PokeDuration = trlPokesMtrx(1,4);
        plxSession(trl).OdorPokeWithdrawTime = trlPokesMtrx(1,2);
        plxSession(trl).MultiOdorPokeLog = false;
        plxSession(trl).OdorPokesDurations = trlPokesMtrx;
        % These check the following
        if ((plxSession(trl).TranspositionDistance==0 && plxSession(trl).Performance==1) ||...
                (plxSession(trl).TranspositionDistance~=0 && plxSession(trl).Performance~=1))
            % Is the trial an InSeq Correct or OutSeq Incorrect? If so,
            % everything's fine, things aren't ambiguous here
        elseif (round(trlPokesMtrx(1,4),4) == round(txtPokeDur(trl),4)) ||...
                round(trlPokesMtrx(1,4),5) == round(txtPokeDur(trl),5)
            % Does the poke withdrawal time match the text file? This
            % happens on InSeq Incorrect and OutSeq Correct. If so,
            % everything's fine, things aren't ambiguous here
        else
            % However, if the poke withdrawal times don't match, that
            % indicates there was an ERROR in assigning the poke duration
            % value in the control code and this trial will require
            % inspection
            plxSession(trl).Errors = [plxSession(trl).Errors; {'ErrType1: No matching poke duration for single poke'}];
            fprintf(outfile, 'ErrType2: No matching poke duration for single poke on trial %i\n', trl);
        end
    %% For trials with multiple pokes
    elseif size(trlPokesMtrx,1)>=2
        plxSession(trl).MultiOdorPokeLog = true;
        %% Step 1: Trial Poke Overview
        % The first check is for instances where everything went correctly,
        if (round(trlPokesMtrx(end,4),4) == round(txtPokeDur(trl),4)) ||...
                (round(trlPokesMtrx(end,4),5) == round(txtPokeDur(trl),5))
            plxSession(trl).PokeDuration = trlPokesMtrx(end,4);
            plxSession(trl).OdorPokeWithdrawTime = trlPokesMtrx(end,2);
            plxSession(trl).OdorPokesDurations = trlPokesMtrx;
            if sum(trlPokesMtrx(1:end-1,5)>= curTxtBuffer)>=1
                plxSession(trl).Errors = [plxSession(trl).Errors; {'ErrType2: Poke buffer error during trial'}];
                fprintf(outfile, 'ErrType2: Poke buffer error; trial continued without reentry during buffer period on trial #%i\n', trl);
            end
        elseif (round(txtPokeDur(trl),4)>=curTargDur && round(trlPokesMtrx(end,4)>=curTargDur)) &&... 
            ((plxSession(trl).TranspositionDistance==0 && plxSession(trl).Performance==1) ||...
            (plxSession(trl).TranspositionDistance~=0 && plxSession(trl).Performance~=1))
            % In these cases the code was working properly, the durations
            % either matched when they should have (withdraw < threshold)
            % or they didn't match when they shouldn't have (withdraw >
            % threshold)
            % To get an accurate withdrawal time, we need to examine the
            % inter-buffer-interval to identify when the rat should have
            % been considered still in the port, and then use that to
            % truncate everything else.
            ipiBufferViolation = find(trlPokesMtrx(:,5) >= curTxtBuffer, 1, 'first');
            if trlPokesMtrx(ipiBufferViolation,4) < curTargDur
                plxSession(trl).Errors = [plxSession(trl).Errors; {'ErrType2: Poke buffer error during trial'}];
                fprintf(outfile, 'ErrType2: Poke buffer error; trial continued without reentry during buffer period on trial #%i\n', trl);
                postThreshPokeOut = find(trlPokesMtrx(:,4) >= curTargDur,1,'first');
                % NOW, shift through the pokes following that one
                % and identify when there was an
                % inter-poke-interval that was greater than the
                % buffer period (i.e. the period where a poke would
                % be considered a continuation). I recognize this is odd,
                % but the buffer is an arbitrary value, this buffer error
                % could reflect them re-positioning within the port, doing
                % things this way gives us some leeway here...
                while trlPokesMtrx(postThreshPokeOut,5) < curTxtBuffer &&...
                        postThreshPokeOut < size(trlPokesMtrx,1)
                    if postThreshPokeOut > size(trlPokesMtrx,1)
                        error('No should go here... wrooooong');
                    end
                    postThreshPokeOut = postThreshPokeOut + 1;
                end
                plxSession(trl).PokeDuration = trlPokesMtrx(postThreshPokeOut,4);
                plxSession(trl).OdorPokeWithdrawTime = trlPokesMtrx(postThreshPokeOut,2);
                plxSession(trl).OdorPokesDurations = trlPokesMtrx(1:postThreshPokeOut,:);
            else
                plxSession(trl).PokeDuration = trlPokesMtrx(ipiBufferViolation,4);
                plxSession(trl).OdorPokeWithdrawTime = trlPokesMtrx(ipiBufferViolation,2);
                plxSession(trl).OdorPokesDurations = trlPokesMtrx(1:ipiBufferViolation,:);
            end
        else
            % For instances where everything did not go perfectly, we need
            % to identify why. Most obvious explanation is that the animal
            % did everything correctly, it just poked during the ITI period
            % and that's being picked up here.
            % The first check to do is to see if there is a difference
            % among the poking events that matches the text file's recorded 
            % poke duration
            %   In instances where poke duration was < threshold these
            %   values should match meaning that 'out' and 'in' will NOT be
            %   empty
            [out, in] = find((round(txtPokeDur(trl),4)==round(trlPokesMtrx(:,2) - trlPokesMtrx(:,1)',4)) |...
                (round(txtPokeDur(trl),5)==round(trlPokesMtrx(:,2) - trlPokesMtrx(:,1)',5)));
            if ~isempty(out) && ~isempty(in)
                % First, we need to validate that the poke IN value is 
                % the start of the trial. If Not, it's an error.
                if in ~= 1
                    plxSession(trl).Errors = [plxSession(trl).Errors; {'ErrType3: Poke initiation was NOT the poke used for timing'}];
                    fprintf(outfile, 'ErrType3: Poke initiation was NOT the poke used for timing, poke #%i was on trial #%i\n', in, trl);
                end                
                plxSession(trl).PokeDuration = trlPokesMtrx(out,4);
                plxSession(trl).OdorPokeWithdrawTime = trlPokesMtrx(out,2);
                plxSession(trl).OdorPokesDurations = trlPokesMtrx(1:out,:);                
                if sum(trlPokesMtrx(1:out-1,5)>= curTxtBuffer)>=1
                    plxSession(trl).Errors = [plxSession(trl).Errors; {'ErrType2: Poke buffer error during trial'}];
                    fprintf(outfile, 'ErrType2: Poke buffer error; trial continued without reentry during buffer period on trial #%i\n', trl);
                end
                if (round(trlPokesMtrx(out,2) - trlPokesMtrx(in,1),4) == round(txtPokeDur(trl),4)) ||...
                        (round(trlPokesMtrx(out,2) - trlPokesMtrx(in,1),5) == round(txtPokeDur(trl),5))
                    % If there is a matching value such that a Poke Out -
                    % Poke In == text poke duration... Perfect! That
                    % accounts for the control code's treatment.
                elseif (round(trlPokesMtrx(out,2) - trlPokesMtrx(in,2),4) == round(txtPokeDur(trl),4)) ||...
                        (round(trlPokesMtrx(out,2) - trlPokesMtrx(in,2),5) == round(txtPokeDur(trl),5))
                    % If there is a matching value BUT it's timed between
                    % poke WITHDRAWAL timepoints.
                    plxSession(trl).Errors = [plxSession(trl).Errors; {'ErrType4: Poke duration was not measured from the poke starting the trial'}];
                    fprintf(outfile, 'ErrType4: Poke duration was measured from poke OUT #%i until poke OUT #%i on trial %i\n', in, out, trl);
                elseif (round(trlPokesMtrx(out,1) - trlPokesMtrx(in,1),4) == round(txtPokeDur(trl),4)) ||...
                        (round(trlPokesMtrx(out,1) - trlPokesMtrx(in,1),5) == round(txtPokeDur(trl),5))
                    % If there is a matching value BUT it's timed between
                    % poke INITIATION timepoints.
                    plxSession(trl).Errors = [plxSession(trl).Errors; {'ErrType5: Poke duration was not measured from the poke starting the trial'}];
                    fprintf(outfile, 'ErrType5: Poke duration was measured from poke IN #%i until poke IN #%i on trial %i\n', in, out, trl);
                elseif (round(trlPokesMtrx(out,1) - trlPokesMtrx(in,2),4) == round(txtPokeDur(trl),4)) ||...
                        (round(trlPokesMtrx(out,1) - trlPokesMtrx(in,2),5) == round(txtPokeDur(trl),5))
                    % If there is a matching value BUT it's timed between
                    % poke WITHDRAWAL and a subsequent poke INITIATION.
                    plxSession(trl).Errors = [plxSession(trl).Errors; {'ErrType6: Poke duration was not measured from the poke starting the trial'}];
                    fprintf(outfile, 'ErrType6: Poke duration was measured from poke OUT #%i until poke IN #%i on trial %i\n', in, out, trl);
                end
            else
                plxSession(trl).Errors = [plxSession(trl).Errors; {'ErrType7: Plexon and Text file poke durations don''t match'}];
                fprintf(outfile, 'ErrType7: Plexon (%.4fs) and Text (%.4fs) poke durations don''t match on trial #%i\n', round(trlPokesMtrx(end,4),5), round(txtPokeDur(trl),5), trl);
                % Here, there is no matching poke duration between the PLX
                % and TXT files, so we need to determine why.
                if round(txtPokeDur(trl),5) < curTargDur
                    % Easiest explanation is that the control code missed
                    % the animal withdrawing and it assigned the standard
                    % value when this contingency was triggered (200ms or
                    % 250ms value depending on the code)
                    % For this, choose the poke withdraw time as the first
                    % time the inter-poke-interval would have violated the
                    % buffer
                    ipiBufferViolation = find(trlPokesMtrx(:,5) >= curTxtBuffer, 1, 'first');
                    plxSession(trl).PokeDuration = trlPokesMtrx(ipiBufferViolation,4);
                    plxSession(trl).OdorPokeWithdrawTime = trlPokesMtrx(ipiBufferViolation,2);
                    plxSession(trl).OdorPokesDurations = trlPokesMtrx(1:ipiBufferViolation,:);
                    plxSession(trl).Errors = [plxSession(trl).Errors; {'ErrType8: Arbitrary withdraw time (likely odor staying on error)'}];
                    fprintf(outfile, 'ErrType8: Arbitrary withdraw time (likely odor staying on error) on trial #%i\n', trl);
                elseif round(txtPokeDur(trl),5) >= curTargDur
                    % The other explanation is that the the trial was held
                    % to, or past, the decision threshold, but he made 
                    % extra pokes after that. Either as a continuation of
                    % his poke, or as ITI pokes.
                    if ((plxSession(trl).TranspositionDistance == 0) &&...
                            (plxSession(trl).Performance == 1))... 
                            ||...
                            ((plxSession(trl).TranspositionDistance ~=0) &&...
                            (plxSession(trl).Performance ~=1))
                        % In this case, use the first poke withdrawal after
                        % the decision threshold as the poke out time
                        postThreshPokeOut = find(trlPokesMtrx(:,4) >= curTargDur,1,'first');
                        % NOW, shift through the pokes following that one
                        % and identify when there was an
                        % inter-poke-interval that was greater than the
                        % buffer period (i.e. the period where a poke would
                        % be considered a continuation)
                        while trlPokesMtrx(postThreshPokeOut,5) < curTxtBuffer &&...
                                postThreshPokeOut < size(trlPokesMtrx,1)
                            if postThreshPokeOut > size(trlPokesMtrx,1)
                                error('No should go here... wrooooong');
                            end
                            postThreshPokeOut = postThreshPokeOut + 1;
                        end                           
                        plxSession(trl).PokeDuration = trlPokesMtrx(postThreshPokeOut,4);
                        plxSession(trl).OdorPokeWithdrawTime = trlPokesMtrx(postThreshPokeOut,2);
                        plxSession(trl).OdorPokesDurations = trlPokesMtrx(1:postThreshPokeOut,:);
                    else
                        error('No code should go here... identify why it''s here')
                    end
                else
                    error('Code no should go here... other cases should be mutually exclusive');
                end                    
            end 
        end
        %% Step 2: Trial Poke Review for potential errors
        % First, review trials where there were additional pokes that
        % occurred PAST the determined withdrawal timepoint. These are
        % easily identified as the matrix stored in the plxSession
        % structure should not match the trlPokesMtrx
        if size(plxSession(trl).OdorPokesDurations,1) ~= size(trlPokesMtrx,1)
            % First instance where this happens is when the buffer was
            % applied incorrectly during the trial
            if plxSession(trl).OdorPokesDurations(end,5) > curTxtBuffer
                % If it was an InSeq Correct, or OutSeq Incorrect trial...
                if ((plxSession(trl).TranspositionDistance == 0) &&...
                        (plxSession(trl).Performance == 1))...
                        ||...
                        ((plxSession(trl).TranspositionDistance ~=0) &&...
                        (plxSession(trl).Performance ~=1))
                    % In this instance, the buffer was run inappropriately
                    % during the trial period, therefore we need to
                    % identify the first poke withdrawal after the decision
                    % threshold
                    postThreshPokeOut = find(trlPokesMtrx(:,4) >= curTargDur,1,'first');
                    % Then, shift through the pokes following that one
                    % and identify when there was an
                    % inter-poke-interval that was greater than the
                    % buffer period (i.e. the period where a poke would
                    % be considered a continuation)
                    while trlPokesMtrx(postThreshPokeOut,5) < curTxtBuffer &&...
                            postThreshPokeOut < size(trlPokesMtrx,1)
                        if postThreshPokeOut > size(trlPokesMtrx,1)
                            error('No should go here... wrooooong');
                        end
                        postThreshPokeOut = postThreshPokeOut + 1;
                    end
                    plxSession(trl).PokeDuration = trlPokesMtrx(postThreshPokeOut,4);
                    plxSession(trl).OdorPokeWithdrawTime = trlPokesMtrx(postThreshPokeOut,2);
                    plxSession(trl).OdorPokesDurations = trlPokesMtrx(1:postThreshPokeOut,:);
                    if postThreshPokeOut ~= find(trlPokesMtrx(:,4) >= curTargDur,1,'first') % This is a little sloppy... but whatever it works                        
                        plxSession(trl).Errors = [plxSession(trl).Errors; {'ErrType10: Buffer error during trial period'}];
                        fprintf(outfile, 'ErrType10: Buffer error during trial period on trial #%i\n', trl);
                    end
                elseif (plxSession(trl).TranspositionDistance ~=0) &&...
                        (plxSession(trl).Performance ==1)...
                        ||...
                        (plxSession(trl).TranspositionDistance ==0) &&...
                        (plxSession(trl).Performance ~=1)
                    % In this case, everything should be ok... but just to
                    % check.
                    if sum((trlPokesMtrx(:,3) > plxSession(trl).OdorPokesDurations(end,4)) &...
                            (trlPokesMtrx(:,3) < curTargDur)) >= 1
                        lateTrialPoke = find((trlPokesMtrx(:,3) > plxSession(trl).OdorPokesDurations(end,4)) &...
                            (trlPokesMtrx(:,3) < curTargDur), 1, 'last');
                        while trlPokesMtrx(lateTrialPoke,5) < curTxtBuffer &&...
                                lateTrialPoke < size(trlPokesMtrx,1)
                            if lateTrialPoke > size(trlPokesMtrx,1)
                                error('No should go here... wrooooong');
                            end
                            lateTrialPoke = lateTrialPoke + 1;
                        end
                        plxSession(trl).PokeDuration = trlPokesMtrx(lateTrialPoke,4);
                        plxSession(trl).OdorPokeWithdrawTime = trlPokesMtrx(lateTrialPoke,2);
                        plxSession(trl).OdorPokesDurations = trlPokesMtrx(1:lateTrialPoke,:);
                        plxSession(trl).Errors = [plxSession(trl).Errors; {'ErrType11: Late trial poke, possibly messy trial'}];
                        fprintf(outfile, 'ErrType11: Late trial poke, possibly messy trial for trial #%i\n', trl);
                    end
                else
                    error('No no spot for code to go.');
                end
            else
%                 if plxSession(trl).OdorPokesDurations(end,4)
                % However, there may just be other buffer issues as well.
                % One example is where the rat triggered the buffer at the
                % end of the trial, but was then back in within the buffer
                % period on the far side of the decision threshold. As this
                % issue is likely him just re-orienting himself, this trial
                % should be considered a potential error.
                % This would have happened on InSeq Incorrect and OutSeq
                % Correct trials
                if (plxSession(trl).TranspositionDistance == 0 && plxSession(trl).Performance ~=1) ||...
                        (plxSession(trl).TranspositionDistance ~=0 && plxSession(trl).Performance ==1)
                    txtTrlEnd = size(plxSession(trl).OdorPokesDurations,1);
                    % Use the same technique to find when the poke SHOULD
                    % have continued until
                    while trlPokesMtrx(txtTrlEnd,5) < curTxtBuffer &&...
                            txtTrlEnd < size(trlPokesMtrx,1)
                        if txtTrlEnd > size(trlPokesMtrx,1)
                            error('No should go here... wrooooong');
                        end
                        txtTrlEnd = txtTrlEnd + 1;
                    end
                    plxSession(trl).PokeDuration = trlPokesMtrx(txtTrlEnd,4);
                    plxSession(trl).OdorPokeWithdrawTime = trlPokesMtrx(txtTrlEnd,2);
                    plxSession(trl).OdorPokesDurations = trlPokesMtrx(1:txtTrlEnd,:);
                    if plxSession(trl).PokeDuration > curTargDur
                        plxSession(trl).Errors = [plxSession(trl).Errors; {'ErrType12: Buffer WOULD have triggered and continued poke past the decision threshold for this trial'}];
                        fprintf(outfile, 'ErrType12: Buffer WOULD have triggered and carried the poke past the decision threshold on trial #%i\n', trl);
                    end                    
                else
                    error('Here too, code should not go.');
                end
            end
                    
        end
    end
    %% Step 3: Identify control code trial attribution errors
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
