function OrganizeStatMatrixByTrial(fileDir)
%% OrganizeStatMatrixByTrial
%   Create a structure with the stat matrix data organized by trial with
%   behavioral data 
%
%   May not be as useful as I initially thought. Doesn't seem to have much
%   more advantage over just extracting timestamps and orienting things
%   with the statMatrix columns... CURRENTLY NOT USED

%% Find directory if none input
if nargin == 0
    fileDir = uigetdir(cd);
    if fileDir == 0
        disp('No directory selected');
        return
    end
    fileDir = [fileDir '\'];
end
fileInfo.OrigDirectory = fileDir;
prompt = {'Pre-Trial Duration (use -0.8 for Irvine, -1 for Boston)', 'Post-Trial Period (use 2 for Irvine, 3 for Boston)'};
dlgTitle = 'Define Trial Period';
numLines = 1;
def = {'-1', '3'};
answer = inputdlg(prompt,dlgTitle,numLines,def);
preWindow = str2double(answer{1});
postWindow = str2double(answer{2});

%% Identify Files
files = dir(fileDir);
fileNames = {files.name};
tetFileLog = cellfun(@(a)~isempty(a), regexp(fileNames,'_T([1-9]*).mat'));
files = fileNames(tetFileLog);
fileInfo.OrigFiles = files;

%% Extract Behavior Data & tsVect
% It's the same in all files so just take it from the first file.
load([fileDir files{1}]);
odorCol = find(cellfun(@(a)~isempty(a), strfind(statMatrixColIDs, 'Odor')), 1, 'first');
behMatrix = statMatrix(:, odorCol:end); %#ok<NODEF>
behMatrixIDs = statMatrixColIDs(odorCol:end);
% Extract tsVect, identfy sample rate & define trial analysis window
tsVect = statMatrix(:,1);
sampleRate = 1/mode(diff(tsVect));
trlWindow = [round(preWindow*sampleRate) round(postWindow*sampleRate)];

%% Create Behavior Vectors for Trial Matrix
% Extract Odor and Position Logical Columns
odorTrlMtx = behMatrix(:,cellfun(@(a)~isempty(a), strfind(behMatrixIDs, 'Odor')));
positionTrlMtx = behMatrix(:,cellfun(@(a)~isempty(a), strfind(behMatrixIDs, 'Position')));
mazePosTrlMtx = behMatrix(:,cellfun(@(a)~isempty(a), strfind(behMatrixIDs, 'MazePosition')));
inSeqVect = behMatrix(:,strcmp(behMatrixIDs, 'InSeqLog'));
perfVect = behMatrix(:,strcmp(behMatrixIDs, 'PerformanceLog'))==1;
frontRewIndx = find(behMatrix(:,strcmp(behMatrixIDs, 'FrontReward')));
rearRewIndx = find(behMatrix(:,strcmp(behMatrixIDs, 'BackReward')));
pokeVect = behMatrix(:,strcmp(behMatrixIDs, 'PokeEvents'));
pokeInIndx = find(pokeVect==1);
pokeOutIndx = find(pokeVect==-1);
trialLog = logical(sum(odorTrlMtx,2));
trialIndx = find(trialLog);
numTrials = sum(trialLog);

trlTS = cell(1,numTrials);
trlStart = cell(1,numTrials);
trlOdors = cell(1,numTrials);
trlPositions = cell(1,numTrials);
trlInSeqLog = cell(1,numTrials);
trlPerf = cell(1,numTrials);
trlFrontRew = cell(1,numTrials);
trlRearRew = cell(1,numTrials);
trlPokeIn = cell(1,numTrials);
trlPokeOut = cell(1,numTrials);
trlPokeDur = cell(1,numTrials);
trlMazePosition = cell(1,numTrials);
trlBlockNumber = cell(1,numTrials);
blockCounter = 0;
for trl = 1:numTrials
    trlOdors{trl} = find(odorTrlMtx(trialIndx(trl),:));
    trlPositions{trl} = find(positionTrlMtx(trialIndx(trl),:)==1,1);
    if trlOdors{trl}==1 && trlPositions{trl}==1
        blockCounter = blockCounter + 1;
    end
    trlBlockNumber{trl} = blockCounter;
    trlInSeqLog{trl} = inSeqVect(trialIndx(trl));
    trlPerf{trl} = perfVect(trialIndx(trl));
    trlPokeIn{trl} = pokeInIndx(find(pokeInIndx<trialIndx(trl), 1,'last'));
    trlPokeOut{trl} = pokeOutIndx(find(pokeOutIndx>trialIndx(trl), 1, 'first'));
    trlPokeDur{trl} = tsVect(trlPokeOut{trl}) - tsVect(trlPokeIn{trl});
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Typically orient relative to Poke IN; make sure to comment in the
    % accompanying fileInfo line when selecting.
    trlStart{trl} = trlPokeIn{trl}; % Use to orient trials to poke in
    fileInfo.TrialStart = 'PokeIn';
%     trlStart{trl} = trialIndx{trl}; % Use to orient trials to odor pres
%     fileInfo.TrialStart = 'OdorPres';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if trl<numTrials
        curFrontRewIndx = frontRewIndx(frontRewIndx>trialIndx(trl) & frontRewIndx<trialIndx(trl+1));
        curRearRewIndx = rearRewIndx(rearRewIndx>trialIndx(trl) & rearRewIndx<trialIndx(trl+1));
    elseif trl==numTrials
        curFrontRewIndx = frontRewIndx(find(frontRewIndx>trialIndx(trl),1,'first'));
        curRearRewIndx = rearRewIndx(find(rearRewIndx>trialIndx(trl),1,'first'));
    end
    trlFrontRew{trl} = tsVect(curFrontRewIndx) - tsVect(trlStart{trl});
    trlRearRew{trl} = tsVect(curRearRewIndx) - tsVect(trlStart{trl});

    curWindow = trlStart{trl} + trlWindow;
    trlMazePosition{trl} = mazePosTrlMtx(curWindow(1):curWindow(2),:);
    trlTS{trl} = tsVect(curWindow(1):curWindow(2)) - tsVect(trlStart{trl});
end

%% Create Data Storage Vectors
trialUnitMatrix = cell(numTrials,length(files));
trialUnitMatrixFR = cell(numTrials,length(files));
trialLFPrawMatrix = cell(numTrials,length(files));
trialLFPhilbMatrix = cell(numTrials,length(files));
unitFRids = [];
unitIDs = [];
tetLFPhilbIDs = [];
tetLFPrawIDs = [];
%% Extract Ensemble Neural Data
for f = 1:length(files)
    load([fileDir files{f}]);
    [tetNameStart, fileTypeStart] = regexp(files{f}, 'T([0-9]*).');
    curTet = files{f}(tetNameStart:fileTypeStart-1);
    % Identify Units
    unitsPosLog = cellfun(@(a)~isempty(a), strfind(statMatrixColIDs, '-U'));
    unitIDs = [unitIDs, tempUnitNames(unitsPosBINlog)]; %#ok<AGROW>
    tempUnitMatrix = statMatrix(:,unitsPosLog);
    
    lfpPosLog = cellfun(@(a)~isempty(a), strfind(statMatrixColIDs, '_LFP'));
    tempLFPnames = statMatrixColIDs(lfpPosLog);
    lfpHilbLog = cellfun(@(a)~isempty(a), strfind(tempLFPnames, '_Hilb'));
    tetLFPhilbIDs = [tetLFPhilbIDs, tempLFPnames(lfpHilbLog)]; %#ok<AGROW>
    lfpRawLog = ~lfpHilbLog;
    tetLFPrawIDs = [tetLFPrawIDs, tempLFPnames(lfpRawLog)]; %#ok<AGROW>
    tempTetLFP = statMatrix(:,lfpPosLog);
    tempTetLFPhilbMatrix = tempTetLFP(:,lfpHilbLog);
    tempTetLFPrawMatrix = tempTetLFP(:,lfpRawLog);
    
    for trl = 1:numTrials
        curWindow = trlStart{trl} + trlWindow;
        
        trialUnitMatrix{trl,f} = tempUnitMatrix(curWindow(1):curWindow(2),:);
        
        trialLFPhilbMatrix{trl,f} = tempTetLFPhilbMatrix(curWindow(1):curWindow(2),:);
        trialLFPrawMatrix{trl,f} = tempTetLFPrawMatrix(curWindow(1):curWindow(2),:);
    end
    fprintf('%s Finished; (%i/%i)\n', curTet, f, length(files));
    clear statMatrix statMatrixColIDs
end

uniEnsemble = cell(1,numTrials);
tetLFPraw = cell(1,numTrials);
tetLFPhilb = cell(1,numTrials);
for trl = 1:numTrials
    uniEnsemble{trl} = cell2mat(trialUnitMatrix(trl,:));
    tetLFPraw{trl} = cell2mat(trialLFPrawMatrix(trl,:));
    tetLFPhilb{trl} = cell2mat(trialLFPhilbMatrix(trl,:));
end
    
fileInfo.UnitIDs = unitIDs;
fileInfo.TetLFPrawIDs = tetLFPrawIDs;
fileInfo.TetLFPhilbIDs = tetLFPhilbIDs;

%% Create statMtxTrlStruct
statMtxTrlStruct = struct('TrialStart', trlStart,'Timestamps', trlTS,...
    'BlockNum', trlBlockNumber, 'Odor', trlOdors, 'Position', trlPositions,...
    'Performance', trlPerf, 'FrontReward', trlFrontRew,...
    'RearReward', trlRearRew, 'InSeq', trlInSeqLog, 'PokeIn', trlPokeIn,...
    'PokeOut', trlPokeOut, 'PokeDur', trlPokeDur, 'MazePosition', trlMazePosition,...
    'EnsembleBIN', uniEnsemble, 'EnsembleFR', uniEnsembleFR,...
    'TetLFPraw', tetLFPraw, 'TetLFPhilb', tetLFPhilb, 'TrialSpectrograms', trialLFPspectMatrix); %#ok<NASGU>

flname = inputdlg('Choose File Name', '1',1,{[fileInfo.OrigDirectory fileInfo.OrigFiles{1}]});
save([cell2mat(flname) '.mat'], 'statMtxTrlStruct', 'fileInfo', '-v7.3');
disp('Trial Structure Saved');