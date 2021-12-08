% IdentifySWR_SM
aniIDs = [{'Barat'},...
    {'Buchanan'},...
    {'Mitt'},...
    {'Stella'},...
    {'SuperChris'}];
dataDir = 'D:\WorkBigDataFiles\CA1 Data\1. WellTrained session\';
aniDirs = cellfun(@(a)sprintf('%s%s\\SWR Tets\\',dataDir, a), aniIDs, 'uniformoutput', 0);
swChans = [{16},...
    {21},...
    {2},...
    {18},...
    {12}];
ripChans = [{21},...
    {18},...
    {18},...
    {14},...
    {15}];
aniInfo = struct('ID', aniIDs, 'Directory', aniDirs, 'SWchan', swChans, 'RIPchan', ripChans,...
    'SWepocs', [], 'SWtrace', [], 'SWpower', [],...
    'RIPepocs', [], 'RIPtrace', [], 'RIPpower', [],...
    'SWRepocs', []);

swThresh = [0 3];
ripThresh = [0 3];
mergeThresh = 15;

odorColors = [44/255, 168/255, 224/255;...
    154/255, 133/255, 122/255;...
    9/255, 161/255, 74/255;...
    128/255, 66/255, 151/255;...
    241/255, 103/255, 36/255];

%% Define Output Data Structure
aniSWRmat = struct('AniID', aniIDs, 'SessionMatrix', [], 'SessionMatrixColIDs', [],...
    'SWRmatrix', [], 'SWRmatrixColIDs', []);
%%
for a = 1:length(aniInfo)
    %% Load Files & Load Behavior Data
    aniFiles = dir(aniInfo(a).Directory);
    behMatFile = aniFiles(cellfun(@(a)~isempty(a), strfind({aniFiles.name}, 'BehaviorMatrix'))).name;
    behav = load([aniInfo(a).Directory behMatFile]);
    behavMat = OrganizeTrialData_SM(behav.behavMatrix, behav.behavMatrixColIDs, [0 0], 'PokeIn');
    sessionMatrix = [[behavMat.PokeInIndex]', [behavMat.PokeOutIndex]', [behavMat.Performance]', [behavMat.Position]', [behavMat.Odor]'];
    aniSWRmat(a).SessionMatrix = sessionMatrix;
    aniSWRmat(a).SessionMatrixColIDs = [{'PokeInIndex'}, {'PokeOutIndex'}, {'Performance'}, {'Trial_Position'}, {'Trial_Odor'}];
    % Extract position and derive velocity 
    posVals = behav.behavMatrix(:,end-1:end);
    locNdxs = find(sum(posVals~=0,2)>=1);
    locVals = behav.behavMatrix(locNdxs,end-1:end)*1.5/100;
    timeVals = behav.behavMatrix(locNdxs,1);
    instV = nan(size(locVals,1),1);
    for v = 2:size(locVals,1)
        instV(v) = sqrt((locVals(v,1)-locVals(v-1,1))^2 - (locVals(v,2)-locVals(v,2))^2)/diff(timeVals(v-1:v));
    end
    smoothInstV = conv(instV, ones(1,20)/20, 'same');
    clear behav behavMat
    
    %% Evaluate Sharpwave channel & events
    swFile = aniFiles(cell2mat(cellfun(@(a)~isempty(a),regexp({aniFiles.name}, ['\w*T' mat2str(aniInfo(a).SWchan) '_\w*SM.mat\>']),'uniformoutput', 0))).name;
    sw = load([aniInfo(a).Directory swFile]);
    [swEpocs, swLFP, swPOW] = SharpwaveDetection(sw.statMatrix(:,2), 1/mode(diff(sw.statMatrix(:,1))), swThresh, mergeThresh);
    swLog = false(size(sw.statMatrix(:,1)));
    for sws = 1:size(swEpocs,1)
        swLog(swEpocs(sws,1):swEpocs(sws,2)) = true;
    end
    swTrace = nan(size(sw.statMatrix(:,1)));
    swTrace(swLog) = sw.statMatrix(swLog,2);
    aniInfo(a).SWepocs = swEpocs;
    aniInfo(a).SWtrace = swLFP;
    aniInfo(a).SWpower = swPOW;
    
    %% Evaluate Ripple channel & events
    ripFile = aniFiles(cell2mat(cellfun(@(a)~isempty(a),regexp({aniFiles.name}, ['\w*T' mat2str(aniInfo(a).RIPchan) '_\w*SM.mat\>']),'uniformoutput', 0))).name;
    rip = load([aniInfo(a).Directory ripFile]);
    [ripEpocs, ripLFP, ripPOW] = RippleDetection(rip.statMatrix(:,2), 1/mode(diff(rip.statMatrix(:,1))), ripThresh, mergeThresh);
    ripLog = false(size(rip.statMatrix(:,1)));
    for rips = 1:size(ripEpocs,1)
        ripLog(ripEpocs(rips,1):ripEpocs(rips,2)) = true;
    end
    ripTrace = nan(size(rip.statMatrix(:,1)));
    ripTrace(ripLog) = rip.statMatrix(ripLog,2);
    aniInfo(a).RIPepocs = ripEpocs;
    aniInfo(a).RIPtrace = ripLFP;
    aniInfo(a).RIPpower = ripPOW;
       
    %% Identify Sharpwave Ripple events
    swrLog = swLog + ripLog;
    swrStart = find(diff(swrLog==2)==1);
    swrEnd = find(diff(swrLog==2)==-1);
    swrWindows = [swrStart nan(length(swrStart),1)];
    for swr = 1:length(swrStart)
        tempSWRstart = swrStart(swr);
        swrWindows(swr,1) = min([swEpocs(tempSWRstart>=swEpocs(:,1) & tempSWRstart<=swEpocs(:,2),1), ripEpocs(tempSWRstart>=ripEpocs(:,1) & tempSWRstart<=ripEpocs(:,2),1), tempSWRstart]);
        
        tempSWRend = swrEnd(find(swrEnd>tempSWRstart,1,'first'));
        swrWindows(swr,2) = max([swEpocs(tempSWRend>=swEpocs(:,1) & tempSWRend<=swEpocs(:,2),2), ripEpocs(tempSWRend>=ripEpocs(:,1) & tempSWRend<=ripEpocs(:,2),2), tempSWRend]);
    end
    for swr = 2:size(swrWindows,1)
        if ~isnan(swrWindows(swr,1)) && ~isnan(swrWindows(swr-1,1))
            if swrWindows(swr,1) - swrWindows(swr-1,2) <= mergeThresh
                swrWindows(swr,1) = swrWindows(swr-1,1);
                swrWindows(swr-1,:) = nan;
            end
        end
    end
    swrWindows(isnan(swrWindows(:,1)),:) = [];
    aniInfo(a).SWRepocs = swrWindows;
    
    newSWtrace = nan(size(sw.statMatrix(:,1)));
    newRIPtrace = nan(size(rip.statMatrix(:,1)));
    for swr = 1:size(swrWindows,1)
        newSWtrace(swrWindows(swr,1):swrWindows(swr,2)) = sw.statMatrix(swrWindows(swr,1):swrWindows(swr,2),2);
        newRIPtrace(swrWindows(swr,1):swrWindows(swr,2)) = rip.statMatrix(swrWindows(swr,1):swrWindows(swr,2),2);
    end
    
    %% Identify Sequence & Non-Sequence SWRs
    intSeqSWRlog = false(size(swrWindows,1),1);
    betSeqsSWRlog = false(size(swrWindows,1),1);
    for swr = 1:size(swrWindows,1)
        preSWRtrlNdx = find(swrWindows(swr,1)>sessionMatrix(1:end-1,2) & swrWindows(swr,1)<sessionMatrix(2:end,1));
        if sessionMatrix(preSWRtrlNdx+1,4)-sessionMatrix(preSWRtrlNdx,4) == 1
            intSeqSWRlog(swr) = true;
        elseif ~isempty(preSWRtrlNdx) && sessionMatrix(preSWRtrlNdx+1,4)==1 && (sessionMatrix(preSWRtrlNdx,4)==max(sessionMatrix(:,4)) || sessionMatrix(preSWRtrlNdx,3)==0)
            betSeqsSWRlog(swr) = true;
        end
    end
%     figure;
%     plot(sessionMatrix(:,1:2)', ones(size(sessionMatrix(:,1:2)')), '-*k');
%     hold on
%     plot(swrWindows(intSeqSWRlog,:)', ones(size(swrWindows(intSeqSWRlog,:)'))-0.1, '-or');
%     plot(swrWindows(betSeqsSWRlog,:)', ones(size(swrWindows(betSeqsSWRlog,:)'))-0.2, '-og');
%     set(gca, 'ylim', [0 1.1]);
    swrMatrix = [swrWindows, intSeqSWRlog, betSeqsSWRlog];
    swrMatrixColIDs = [{'SWRstartIndex'}, {'SWRendIndex'}, {'WithinSeqSWRlog'}, {'BetweenSeqSWRlog'}];
    aniSWRmat(a).SWRmatrix = swrMatrix;
    aniSWRmat(a).SWRmatrixColIDs = swrMatrixColIDs;
    
    
    save([dataDir aniIDs{a} '\' aniIDs{a} '_SWRs.mat'], 'swrMatrix', 'swrMatrixColIDs');

    %% Plot the session SWR traces
    figure; 
    sp1 = subplot(2,1,1);
    plot(swTrace);
    hold on; 
    plot(ripTrace);
    title(aniInfo(a).ID);
      
    sp2 = subplot(2,1,2);
    plot(newSWtrace);
    hold on; 
    plot(newRIPtrace);
    
    linkaxes([sp1 sp2], 'xy');
    
    annotation(gcf,'textbox', [0 0.01 1 0.05],'String', ['SW File:' swFile],...
        'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
    annotation(gcf,'textbox', [0 0.05 1 0.05],'String', ['RIP File: ' ripFile],...
        'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
end

%% %%%%%%%%%
%  Functions
%  %%%%%%%%%

%%
function [ripEpocs, lfpFilt, zPow] = RippleDetection(lfp, fs, powThresh, mergeThresh)
%% RippleDetection
% Code for detection of ripples on tetrode presumed to be in pyramidal
%   layer of CA1
% Inputs:
%   lfp: Nx1 input of the raw LFP 
%   powThresh: 1x2 input with threshold parameters for detection of peak
%       event. Initial value is event boundaries, second is peak required
%   mergeThresh: threshold for merging two sharpwave events if they happen
%       in close temporal proximity%
%%
if nargin == 0
    [file,path] = uigetfile('*.mat');
    load([path file]); %#ok<LOAD>
    ripCol = cellfun(@(a)~isempty(a), regexp(statMatrixColIDs, 'LFP_Ripple$'));
    lfp = statMatrix(:,ripCol); %#ok<USENS>
    clear statMatrix statMatrixColIDs file path
    powThresh = [0 4];
    mergeThresh = 15;
end
%%
lfpFilt = SimpleFilt(lfp, fs, 150, 250);
zPow = zscore(abs(hilbert(lfpFilt)));

%% Identify Potential Sharpwave Periods
abvThresh1 = zPow > powThresh(1);
epocStart = find(diff(abvThresh1)==1);
epocEnd = find(diff(abvThresh1)==-1);
epocWindows = [epocStart, nan(size(epocStart))];
for epoc = 1:size(epocStart,1)
    tempEpocEnd = epocEnd(find(epocEnd>epocStart(epoc), 1, 'first'));
    if sum(zPow(epocStart(epoc):tempEpocEnd)>=powThresh(2))>=1
        epocWindows(epoc,2) = tempEpocEnd;
    else
        epocWindows(epoc,1) = nan;
    end
end
epocWindows(isnan(epocWindows(:,1)),:) = [];
% Merge short latency sharpwaves
for epoc = 2:size(epocWindows,1)
    if ~isnan(epocWindows(epoc,1)) && ~isnan(epocWindows(epoc-1,1))
        if epocWindows(epoc,1) - epocWindows(epoc-1,2) <= mergeThresh
            epocWindows(epoc,1) = epocWindows(epoc-1,1);
            epocWindows(epoc-1,:) = nan;
        end
    end
end
epocWindows(isnan(epocWindows(:,1)),:) = [];
%%
ripEpocs = epocWindows;
end

function [swEpocs, hpf, zPow] = SharpwaveDetection(lfp, fs, powThresh, mergeThresh)
%% SharpwaveDetection
% Code for detection of sharpwaves on tetrodes presumed to be in the 
%   stratum radiatum of CA1
% Inputs:
%   lfp: Nx1 input of the raw LFP 
%   powThresh: 1x2 input with threshold parameters for detection of peak
%       event. Initial value is event boundaries, second is peak required
%   mergeThresh: threshold for merging two sharpwave events if they happen
%       in close temporal proximity
%%
if nargin == 0
    [file,path] = uigetfile('*.mat');
    load([path file]); %#ok<LOAD>
    lfp = statMatrix(:,2); %#ok<USENS>
    fs = 1/mode(diff(statMatrix(:,1)));
    clear statMatrix statMatrixColIDs file path
    powThresh = [0 2];
    mergeThresh = 15;
end
%% Filter & Z-Score Envelope
hpf = highpass(lfp, 4, fs);
zPow = zscore(abs(hilbert(hpf)));

%% Identify Potential Sharpwave Periods
abvThresh1 = zPow > powThresh(1);
epocStart = find(diff(abvThresh1)==1);
epocEnd = find(diff(abvThresh1)==-1);
epocWindows = [epocStart, nan(size(epocStart))];
for epoc = 1:size(epocStart,1)
    tempEpocEnd = epocEnd(find(epocEnd>epocStart(epoc), 1, 'first'));
    if sum(zPow(epocStart(epoc):tempEpocEnd)>=powThresh(2))>=1
        epocWindows(epoc,2) = tempEpocEnd;
    else
        epocWindows(epoc,1) = nan;
    end
end
epocWindows(isnan(epocWindows(:,1)),:) = [];
% Merge short latency sharpwaves
for epoc = 2:size(epocWindows,1)
    if ~isnan(epocWindows(epoc,1)) && ~isnan(epocWindows(epoc-1,1))
        if epocWindows(epoc,1) - epocWindows(epoc-1,2) <= mergeThresh
            epocWindows(epoc,1) = epocWindows(epoc-1,1);
            epocWindows(epoc-1,:) = nan;
        end
    end
end
epocWindows(isnan(epocWindows(:,1)),:) = [];
%%
swEpocs = epocWindows;
%%
end

%% Simple Filtering
function lfpFilt = SimpleFilt(trace, samp, low, high)
Wn_FRange = [low/(samp/2) high/(samp/2)]; % normalized by the nyquist frequency
[bFRange, aFRange] = butter(3,Wn_FRange);
lfpFilt = filtfilt(bFRange,aFRange,trace);
end