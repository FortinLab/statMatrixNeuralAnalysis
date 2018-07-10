function [statMatrixColIDs, statMatrix] = CreateStatMatrixFromOE(recStart,processorNum,channelNum,fileDir,ssnDir,varargin)
%% Updated 7/10/2018: Added 'PokeInEvents' Column and information for additionional frequency bands. Updated for new EventsTS script(formerly BehavTS)
%
% This is a script for creating a statMatrix for a single .continuous file
% Example input for running statMatrix script for Channel 1 of a session:
% [statMatrixColIDs,statMatrix] = CreateStatMatrixFromOE_new(recStart, 100 ,1, fileDir, ssnDir,'processor', 'ADC', 'CH','session','allfreq')
% ProcessorNum = 100 (FPGA or Raw); processorNum = 102 (Bandpass)
% fileDir = the folder directory of .continuous, .events, and/or .spike files
% For recStart, open the messages.events file (using Notepad or Notepad++) in target file directory to find starting time of recording
% ssnDir = the folder directory and the file name of your session data

%% Argument String Detection
% Determines what data to extract from data sets

if sum(cell2mat(cellfun(@(a)(strcmp(a,'ADC') | strcmp(a, 'adc')), varargin, 'uniformoutput', 0)))>=1
    [Events_TS,~,~,~] = EventsTS(fileDir,processorNum,recStart,'downsample');
end

if sum(cell2mat(cellfun(@(a)(strcmp(a,'CH') | strcmp(a, 'ch')), varargin, 'uniformoutput', 0)))>=1
    filename = sprintf('%d_CH%d.continuous', processorNum,channelNum);
    [LFP_Raw,~,info] = load_open_ephys_data_faster([fileDir '\' filename]);
    sampleRate = info.header.sampleRate;
    TS_Raw = recStart:1/sampleRate:recStart+(length(LFP_Raw)/sampleRate)-(1/sampleRate);
end

if sum(cell2mat(cellfun(@(a)(strcmp(a, 'session') | strcmp(a, 'Session')), varargin, 'uniformoutput', 0)))>=1
    load(ssnDir,'ssnData')
end

% Decision call for band pass filtering up to a certain frequency band
% If I call 'theta', it only filters for theta.  If I call 'beta', it
% filters up to beta aka theta, low beta, and beta.

if sum(cell2mat(cellfun(@(a)(strcmp(a, 'allFreq') | strcmp(a, 'allfreq') | strcmp (a,'Allfreq')...
        | strcmp(a,'AllFreq')), varargin, 'uniformoutput', 0)))>=1
    numfreqBand = 6; % For bandpass filtering all frequency bands
elseif sum(cell2mat(cellfun(@(a)(strcmp(a, 'Theta') | strcmp(a, 'theta')), varargin, 'uniformoutput', 0)))>=1
    numfreqBand = 1;
elseif sum(cell2mat(cellfun(@(a)(strcmp(a, 'Low beta') | strcmp(a, 'Low Beta')), varargin, 'uniformoutput', 0)))>=1
    numfreqBand = 2;
elseif sum(cell2mat(cellfun(@(a)(strcmp(a, 'Beta') | strcmp(a, 'beta')), varargin, 'uniformoutput', 0)))>=1
    numfreqBand = 3;
elseif sum(cell2mat(cellfun(@(a)(strcmp(a, 'Low gamma') | strcmp(a, 'Low Gamma')), varargin, 'uniformoutput', 0)))>=1
    numfreqBand = 4;
elseif sum(cell2mat(cellfun(@(a)(strcmp(a, 'High gamma') | strcmp(a, 'High Gamma')), varargin, 'uniformoutput', 0)))>=1
    numfreqBand = 5;
end

%% Bandpass Filter and Hilbert Transform of LFP

% Lowpass Filtered 0-300Hz and Transformed

order = 3;
fcut  = 300;
NyquistFreq = sampleRate/2;
[b,l] = butter(order,fcut/(NyquistFreq),'low');
LFP_Raw = filter(b,l,LFP_Raw);

LFPHilb_Raw = atan2(imag(hilbert(LFP_Raw)), LFP_Raw); % Corresponding Hilbert Transform

% Check if signal needs to be downsampled before the Hilbert transforms
if sampleRate > 1000
    n_samples = sampleRate/1000;
    LFP_Raw = downsample(LFP_Raw,n_samples);
    LFPHilb_Raw = downsample(LFPHilb_Raw,n_samples);
    TS_Raw = downsample(TS_Raw,n_samples);
end

% Bandpass filtering for all frequencies and their Hilbert values
for freqBand = 1:numfreqBand
    if freqBand == 1
        MinFreq = 4;
        MaxFreq = 12;
        wSize = 1.25*(MaxFreq+MinFreq)/2;
        [LFPHilb_Theta, LFP_Theta] = PhaseFreqDetectAbbr(LFP_Raw,TS_Raw,MinFreq,MaxFreq);
    end
    if freqBand == 2
        MinFreq = 13;
        MaxFreq = 19;
        wSize = 1.25*(MaxFreq+MinFreq)/2;
        [LFPHilb_LowBeta, LFP_LowBeta] = PhaseFreqDetectAbbr(LFP_Raw,TS_Raw,MinFreq,MaxFreq);
    end
    if freqBand == 3
        MinFreq = 20;
        MaxFreq = 40;
        wSize = 1.25*(MaxFreq+MinFreq)/2;
        [LFPHilb_Beta, LFP_Beta] = PhaseFreqDetectAbbr(LFP_Raw,TS_Raw,MinFreq,MaxFreq);
    end
    if freqBand == 4
        MinFreq = 41;
        MaxFreq = 59;
        wSize = 1.25*(MaxFreq+MinFreq)/2;
        [LFPHilb_LowGamma, LFP_LowGamma] = PhaseFreqDetectAbbr(LFP_Raw,TS_Raw,MinFreq,MaxFreq);
    end
    if freqBand == 5
        MinFreq = 61;
        MaxFreq = 80;
        wSize = 1.25*(MaxFreq+MinFreq)/2;
        [LFPHilb_HighGamma, LFP_HighGamma] = PhaseFreqDetectAbbr(LFP_Raw,TS_Raw,MinFreq,MaxFreq);
    end
    if freqBand == 5
        MinFreq = 150;
        MaxFreq = 250;
        wSize = 1.25*(MaxFreq+MinFreq)/2;
        [LFPHilb_Ripple, LFP_Ripple] = PhaseFreqDetectAbbr(LFP_Raw,TS_Raw,MinFreq,MaxFreq);
    end
end



%% Organize LFP Columns
% Organizes all LFP related vectors into columns

TS_Col = transpose(TS_Raw);

LFP_Mat = [LFP_Raw,LFPHilb_Raw,LFP_Theta,LFPHilb_Theta, LFP_LowBeta, LFPHilb_LowBeta, LFP_Beta, LFPHilb_Beta,...
            LFP_LowGamma, LFPHilb_LowGamma, LFP_HighGamma, LFPHilb_HighGamma, LFP_Ripple, LFPHilb_Ripple];

LFPDataHeaders = ["LFP_Raw","LFP_Raw_Hilbert","LFP_Theta","LFP_Theta_Hilbert","LFP_LowBeta","LFP_LowBeta_Hilbert",...
    "LFP_Beta","LFP_Beta_Hilbert","LFP_LowGamma","LFP_LowGamma_Hilbert","LFP_HighGamma","LFP_HighGamma_Hilbert",...
    "LFP_Ripple","LFP_Ripple_Hilbert"];


%% Create and Organize Behavioural Columns

% Initializes Behavioral Column Headers for Label References
seqLength = ssnData(1).Settings.SequenceLength;

for seq = 1:seqLength
    behDataHeaders{seq} = ['Odor' num2str(seq)];
    behDataHeaders{seq+seqLength} = ['Position' num2str(seq)];
end

behDataHeaders = [behDataHeaders, "InSeqLog", "Performance", "PokeInEvents"];

% Odor onset data from ADC .continuous files
Odor1_Col = transpose(Events_TS.OnPeakTS_1);
Odor2_Col = transpose(Events_TS.OnPeakTS_2);
Odor3_Col = transpose(Events_TS.OnPeakTS_3);
Odor4_Col = transpose(Events_TS.OnPeakTS_4);

% Extracts indices for odor onsets
Odor1_Index = find(Events_TS.OnPeakTS_1);
Odor2_Index = find(Events_TS.OnPeakTS_2);
Odor3_Index = find(Events_TS.OnPeakTS_3);
Odor4_Index = find(Events_TS.OnPeakTS_4);

% Combines odor events indices for use in later columns
AllOdorOnset = Events_TS.OnPeakTS_1 + Events_TS.OnPeakTS_2 +...
    Events_TS.OnPeakTS_3 + Events_TS.OnPeakTS_4; 
AllOdorOnset_Index = find(AllOdorOnset==1);

% Extracts Poke-in events data

[~,PokeInEventsIndices] = intersect(TS_Raw,OnTS.OnPeakTS_5);
PokeInEvents_Col = zeroes(length(TS_Col),1);
PokeInEvents_Col(PokeInEventsIndices) = 1;

% Organizes Performance data from session file
Performance_Session = [ssnData.Performance];
for PerformanceIndex = 1:size(AllOdorOnset_Index,2)
    if Performance_Session(PerformanceIndex) == AllOdorOnset(AllOdorOnset_Index(PerformanceIndex))
        AllOdorOnset(AllOdorOnset_Index(PerformanceIndex)) = 1;
    else
        AllOdorOnset(AllOdorOnset_Index(PerformanceIndex)) = -1;
    end
end
Performance_Col = transpose(AllOdorOnset);

% Organizes In-sequence trials data from session file
SessionLog = [ssnData.TranspositionDistance];
InSeqLog = zeros(1,size(Events_TS.OnPeakTS_1,2));
for SessionLogIndex = 1:size(SessionLog,2)
    if SessionLog(SessionLogIndex) == 0
        InSeqLog(AllOdorOnset_Index(SessionLogIndex)) = 1;
    else
        InSeqLog(AllOdorOnset_Index(SessionLogIndex)) = -1;
    end
end
InSeqLog_Col = transpose(InSeqLog);

% Organizes Position data from Odor Indices
Position1 = histcounts(TS_Raw(Odor1_Index),[0 TS_Raw]);
Position2 = histcounts(TS_Raw(Odor2_Index),[0 TS_Raw]);
Position3 = histcounts(TS_Raw(Odor3_Index),[0 TS_Raw]);
Position4 = histcounts(TS_Raw(Odor4_Index),[0 TS_Raw]);

Position1_Col = transpose(Position1);
Position2_Col = transpose(Position2);
Position3_Col = transpose(Position3);
Position4_Col = transpose(Position4);

% Creates Behavrioral Data Matrix 
Behav_Mat = [Odor1_Col,Odor2_Col,Odor3_Col,Odor4_Col,Position1_Col,Position2_Col,...
    Position3_Col,Position4_Col,InSeqLog_Col,Performance_Col,PokeInEvents_Col];

%% Arrange all matrices/columns into full statMatrix

statMatrix = [TS_Col LFP_Mat Behav_Mat];
statMatrixColIDs = ['TimeBins', LFPDataHeaders, behDataHeaders];

end



