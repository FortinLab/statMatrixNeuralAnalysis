function [ssnData] = FillPLXeventTSsToSsnData(plxFile, matFile)
%% FillPLXeventTSsToSsnData
%   Alternative approach to examining timestamps/trial values etc. This one
%   uses the parameters stored in the ssnData .mat file to guide the
%   operations of the timestamp extraction...
%
%   06/28/2018 - Created by GE
%% Identify files
clc
if ~(nargin == 2)
    if nargin == 0
        [fileName, path] = uigetfile('.plx','Identify .PLX File');
        if fileName == 0
            disp('No file selected, analysis cancelled')
            return
        end
        plxFile = [path fileName];
    end
    [path, filename] = fileparts(plxFile);
    path = [path '\'];
    flContents = dir(path);
    fileNames = {flContents.name};
    matFileLog = cellfun(@(a)~isempty(a), regexp(fileNames, [filename '_([0-9]*)-([A-Z | a-z]*)-([0-9]*).mat']));
    matFile = [path fileNames{matFileLog}];    
end

%% Load Data
% Load the ssnData file first
load(matFile);
if isfield(ssnData(1).Settings, 'ShortPokeBufferDur')
    bufferPeriod = ssnData(1).Settings.ShortPokeBuffer;
    bufferDuration = ssnData(1).Settings.ShortPokeBufferDur;
else
    error('This data was collected with an old version of the code that didn''t allow declaration of buffer duration, please use an older plxAnalysis script for it');
end
% Now extract event data from the .plx file
[numChans, chanNames] = plx_event_names(plxFile);
plxStruct = struct('Channel', [], 'n', [], 'ts', [], 'sv', []);
plxSsn = cell(numChans,1);
for chan = 1:numChans
    curChan = chanNames(chan,:);
    intVals = double(curChan);
    valLim = find(~(intVals==0), 1, 'last');
    plxStruct(chan).Channel = curChan(1:valLim);
    [plxStruct(chan).n, plxStruct(chan).ts, plxStruct(chan).sv] = plx_event_ts(plxFile, curChan);
    plxSsn{chan} = [plxStruct(chan).ts, ones(size(plxStruct(chan).ts)).*chan];
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
odorPresTime = cell(1,10); % ASSUMES 10 ODORS; CHANGE IF THIS ASSUMPTION IS NO LONGER CORRECT, i.e. if there are >10

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

odorPresSsn = odorPresTime;
for opt = 1:length(odorPresSsn)
    odorPresSsn{opt} = [ones(length(odorPresSsn{opt}),1).*opt odorPresSsn{opt}];
end
odorPresSsn(cellfun(@(a)isempty(a), odorPresSsn)) = [];
odorPresSsn = sortrows(cell2mat(odorPresSsn'), 2);

if ~sum(odorPresSsn(:,1)'-[ssnData.Odor])==0
    fprintf('PLX file = %s\n', plxFile);
    fprintf('MAT file = %s\n', matFile);
    error('Odor order in .plx and .mat files don''t align, ensure you are working with the correct files')
end

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
interPokeInterval(1) = []; % Remove first inter-poke interval because it refers to period between first and second poke, not period from start to end of first poke
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
% Something is wrong with the way these flags are inserted into the
% recordings, for now it's not important since we're combining with the
% ssnData value anyway it's fine to ignore these events channels and use
% the ssnData structure for this info.

% corrTrialChanNum = strcmp(channels, 'Correct Poke') | strcmp(channels, 'Keyboard7');
% correctTrialsTimes = plxStruct(corrTrialChanNum).ts;
% 
% inCorrTrialChanNum = strcmp(channels, 'Incorrect Poke') | strcmp(channels, 'Keyboard8');
% incorrectTrialsTimes = plxStruct(inCorrTrialChanNum).ts;
% 
% if term
%     correctTrialsTimes(correctTrialsTimes>terminateTime) = [];
%     incorrectTrialsTimes(incorrectTrialsTimes>terminateTime) = [];
% end
% 
% trialOutcomes = sort([correctTrialsTimes;incorrectTrialsTimes]);
% if sum(pokeInitiationTimes(1) > trialOutcomes) >= 1 % If a trial outcome is registered before a poke is registered that means there was a poke already in the buffer and matlab saw that poke in but not the accompanying poke out. This happens when the nose port is cleaned between sessions and the buffer is not cleared before running the next session.
%     miscreant = trialOutcomes(pokeInitiationTimes(1) > trialOutcomes);
%     correctTrialsTimes(correctTrialsTimes==miscreant) = [];
%     incorrectTrialsTimes(incorrectTrialsTimes==miscreant) = [];
% end
% perfSsn = sortrows([correctTrialsTimes ones(size(correctTrialsTimes)); incorrectTrialsTimes ones(size(incorrectTrialsTimes))]);

%% Trial End Times
% This flag is inserted at the start and end of every trial... it should
% probably be called Trial Boundaries come to think of it... too late....
trialEndChanNum = strcmp(channels, 'Trial End');
trialEndTimes = plxStruct(trialEndChanNum).ts;

if term
    trialEndTimes(end) = [];
end

if ~(length(trialEndTimes)==length(ssnData)*2)
    fprintf('PLX file = %s\n', plxFile);
    fprintf('MAT file = %s\n', matFile);
    error('More Trial End Time events than trials in ssnData, check files');
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

%% Extract Session Timestamps
trl = 1;

