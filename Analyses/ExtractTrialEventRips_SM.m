function [trialRipStruct] = ExtractTrialEventRips_SM(rips, trialWin)
% ExtractTrialEventRips_SM
%   Extracts and organizes ripple data relative to the trial periods
%
%   Inputs:
%       - rips: Ripple data structure created by RippleDetection_SM
%       - trialWin: Periods used for delineation of trial event related
%           SWR.
%               - First value is period prior to poke initiation considered
%                   as the "pre-trial" period
%               - Second value is the period after poke withdrawal
%                   considered as the "post-trial" period

%% Extract Near-Trial Epocs
trialPokeTimes = rips.TrialInfo.TrialPokes;
trialRips = cell(size(trialPokeTimes,1),3);
trialRipLat = cell(size(trialPokeTimes,1),3);
trialRipsDur = cell(size(trialPokeTimes,1),3);
trialRipsSync = cell(size(trialPokeTimes,1),3);
for trl = 1:size(trialPokeTimes,1)
    preTrlLog = rips.Ripples.Events(:,1)>(trialPokeTimes(trl,1)-trialWin(1)) & rips.Ripples.Events(:,1)<trialPokeTimes(trl,1);
    trialRips{trl,1} = rips.Ripples.Events(preTrlLog,:);
    trialRipLat{trl,1} = rips.Ripples.Events(preTrlLog,1) - trialPokeTimes(trl,1);
    trialRipsDur{trl,1} = rips.Ripples.Duration(preTrlLog,:);
    trialRipsSync{trl,1} = rips.Ripples.Synchrony(preTrlLog,:);
    
    trlLog = rips.Ripples.Events(:,1)>trialPokeTimes(trl,1) & rips.Ripples.Events(:,1)<trialPokeTimes(trl,2);
    trialRips{trl,2} = rips.Ripples.Events(trlLog,:);
    trialRipLat{trl,2} = rips.Ripples.Events(trlLog,1) - trialPokeTimes(trl,1);
    trialRipsDur{trl,2} = rips.Ripples.Duration(trlLog,:);
    trialRipsSync{trl,2} = rips.Ripples.Synchrony(trlLog,:);
    
    pstTrlLog = rips.Ripples.Events(:,1)>trialPokeTimes(trl,2) & (rips.Ripples.Events(:,1)<trialPokeTimes(trl,2)+trialWin(2));
    trialRips{trl,3} = rips.Ripples.Events(pstTrlLog,:);
    trialRipLat{trl,3} = rips.Ripples.Events(pstTrlLog,1) - trialPokeTimes(trl,2);
    trialRipsDur{trl,3} = rips.Ripples.Duration(pstTrlLog,:);
    trialRipsSync{trl,3} = rips.Ripples.Synchrony(pstTrlLog,:);    
end

trialRipStruct = struct('Events', {trialRips}, 'Latency', {trialRipLat}, 'Duration', {trialRipsDur}, 'Synchrony', {trialRipsSync});