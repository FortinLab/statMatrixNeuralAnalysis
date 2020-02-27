function [rips] = RippleDetection_SM(envProc, powThresh, durThresh, durThreshMrg, syncThresh, syncWin, smoothWin)
%% RippleDetection_SM
%   Inputs:
%       - envProc: Enveloping procedure
%               - 'RMS' : Use the root mean squared approach to enveloping
%               - 'HILB' : Use hilbert transform envelope
%       - powThresh: Thresholds used for defining the ripple.
%               - First value is the threshold for defining a ripple epoc.
%               - Second value is the threshold needed to for an epoc to be
%                   considered a ripple.
%       - durThresh: Threshold used to select ripples of only a particular
%           duration. **CURRENTLY NOT IMPLEMENTED**
%       - durThreshMrg: Threshold of the inter-ripple-interval used to
%           merge together potential ripples that happen close in time
%       - syncThresh: Synchrony threshold used to select ripples of a
%           particular coherence. **CURRENTLY NOT IMPLEMENTED**
%       - syncWin: Window size used to compute sliding window coherence
%           across the session. **CURRENTLY NOT IMPLEMENTED**
%       - smoothWin: Window size used for gaussian smoothing.


%#ok<*IDISVAR,*NODEF,*USENS,*NASGU,*COLND>
%%
if nargin == 0
    envProc = 'RMS';            % Enable for RMS determination of envelope
    % envProc = 'HILB';           % Enable for abs(Hilbert) determinination of envelope
    powThresh = [0 4];
    durThresh = 15;             % Duration Threshold
    durThreshMrg = 15;
    syncThresh = 0;
    syncWin = 10;
    smoothWin = 21;
end
%%
%% Define Analysis Features
w = gausswin(smoothWin);
w = w/sum(w);
%% Load Data
origDir = cd;
[fileDir] = uigetdir(origDir);
if fileDir==0
    disp('Analysis Cancelled')
    return
else
    cd(fileDir)
end
dirContents = dir(fileDir);
fileNames = {dirContents.name};
tetFileLog = cellfun(@(a)~isempty(a), regexp(fileNames, '_T([0-9]*)_SM.mat'));
tetFiles = fileNames(tetFileLog)';
behavFileLog = cellfun(@(a)~isempty(a), regexp(fileNames, '_BehaviorMatrix'));
ensembleFileLog = cellfun(@(a)~isempty(a), regexp(fileNames, '_EnsembleMatrix'));
load(fileNames{behavFileLog});
load(fileNames{ensembleFileLog});
behavMatrixTrialStruct = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [0 0], 'PokeIn');

%% Extract Raw Values & Compute RMS Power
ripBPF = nan(size(behavMatrix,1),size(tetFiles,1));
ripVolts = nan(size(behavMatrix,1),size(tetFiles,1));
ripRMS = nan(size(behavMatrix,1),size(tetFiles,1));
ripHilb = nan(size(behavMatrix,1),size(tetFiles,1));
ripTetIDs = cell(size(tetFiles));
for fl = 1:length(tetFiles)
    load(tetFiles{fl})
    samp = mode(diff(statMatrix(:,1)));
    wIndx = round((1/200)*(1/samp));
    fprintf('%s......', tetFiles{fl});
    ripCol = find(cellfun(@(a)~isempty(a), regexp(statMatrixColIDs, 'LFP_Ripple$')));
    ripVolts(:,fl) = statMatrix(:,2);
    if strcmp(envProc, 'RMS')
        ripRMS(:,fl) = conv(sqrt(conv(statMatrix(:,ripCol).^2, ones(wIndx,1)/wIndx, 'same')), w, 'same');
    elseif strcmp(envProc, 'HILB')
        ripRMS(:,fl) = conv(abs(hilbert(statMatrix(:,ripCol))),w,'same');
    end
    ripBPF(:,fl) = statMatrix(:,ripCol);
    ripHilb(:,fl) = statMatrix(:,ripCol+1);
    ripTetIDs{fl} = statMatrixColIDs{ripCol};
    fprintf('done\n');
end

%% Calculate Thresholds
% Aggregate Power
aggPower = mean(ripRMS,2);                  % Mean envelope
% aggPower = zscore(mean(ripRMS,2));          % Z-Score Mean envelope
% aggPower = mean(zscore(ripRMS),2);          % Mean Z-Score envelope

% Threshold Based on Mean +/- STD Aggregate Power
rmsThresh1 = (mean(aggPower) + (powThresh(1)*std(aggPower)));
rmsThresh2 = (mean(aggPower) + (powThresh(2)*std(aggPower)));

%% Identify Ripples
% Define putative ripple periods
abvThresh1 = aggPower>rmsThresh1;
epocWindows = [find(diff(abvThresh1)==1), find(diff(abvThresh1)==-1)];

% Apply secondary (peak power) threshold
abvThresh2 = aggPower>rmsThresh2;
dualThreshLog = false(size(epocWindows,1),1);
for epoc = 1:size(epocWindows,1)
    if sum(abvThresh2(epocWindows(epoc,1):epocWindows(epoc,2))) >=1
        dualThreshLog(epoc) = true;
    end
end
epocWindows(~dualThreshLog,:) = [];                                         % Comment out if not using dual thresholds

% Merge short latency ripples
interEpocInterval = epocWindows(2:end,1) - epocWindows(1:end-1,2);
slrLog = interEpocInterval<durThreshMrg;
shortLatRipls = find(slrLog);
mrgdNdxs = false(size(interEpocInterval));
for slr = 1:length(shortLatRipls)
    nxtNdx = shortLatRipls(slr)+find(slrLog(shortLatRipls(slr)+1:end)==0,1,'first');
    epocWindows(shortLatRipls(slr),2) = epocWindows(nxtNdx,2);
    mrgdNdxs(nxtNdx) = true;
end
epocWindows(mrgdNdxs,:) = [];

%% Quantify & Extract Ripples
% Determine the duration/timing of each event
epocDur = diff(epocWindows,1,2);

% Determine the Synchrony of each event
epocSync = nan(size(epocWindows,1),1);
epocSyncConfMtx = cell(size(epocWindows,1),1);
for epoc = 1:size(epocWindows,1)
    tempConfMtx = nan(size(ripHilb,2));
    for e1 = 1:size(ripHilb,2)
        for e2 = 1:size(ripHilb,2)
            [tempConfMtx(e1,e2),~] = circ_corrcc(ripHilb(epocWindows(epoc,1):epocWindows(epoc,2),e1),...
                ripHilb(epocWindows(epoc,1):epocWindows(epoc,2),e2));
        end
    end
    epocSyncConfMtx{epoc} = tempConfMtx;
    epocSync(epoc) = mean(tempConfMtx(triu(true(size(ripHilb,2)))));
end

%% Organize Spiking Data Based on Total # Spikes
spkMtx = ensembleMatrix(:,2:end);
sortedSpkCountsAndIndices = sortrows([sum(spkMtx); 1:size(spkMtx,2)]');
spkMtxSorted = spkMtx(:,sortedSpkCountsAndIndices(:,2));


%% Organize Data Output
rips = struct(...
    'TimeStamps', statMatrix(:,1),...
    'Ripples', struct('Events', epocWindows, 'Duration', epocDur,...
        'Synchrony', epocSync),...
    'SessionData', struct('RawLFP', ripVolts, 'RipBPF', ripBPF,...
        'RipEnv', ripRMS, 'RipPhase', ripHilb, 'TetIDs', {ripTetIDs},...
        'Spikes', spkMtxSorted),...
    'TrialInfo', struct('Perf', [behavMatrixTrialStruct.Performance]==1,...
        'TransDist', [behavMatrixTrialStruct.TranspositionDistance],...
        'OdorVect', [behavMatrixTrialStruct.Odor],...
        'PositionVect', [behavMatrixTrialStruct.Position],...
        'TrialPokes', [[behavMatrixTrialStruct.PokeInIndex]', [behavMatrixTrialStruct.PokeOutIndex]'],...
        'TrialRewards', [behavMatrixTrialStruct.RewardIndex]),...
    'FileInfo', struct('Directory', cd, 'EnvelopeProcedure', envProc,...
        'PowerThreshold', powThresh, 'DurationThreshold', durThresh,...
        'DurationMergeThreshold', durThreshMrg,...
        'SynchronyThreshold', syncThresh, 'SynchronyWindow', syncWin,...
        'GaussianDuration', smoothWin));
end