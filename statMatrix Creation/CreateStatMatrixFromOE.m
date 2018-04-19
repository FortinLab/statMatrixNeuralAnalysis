function [statMatrixColIDs, statMatrix, AnalysisColumns] = CreateStatMatrixFromOE(recStart,processorNum,channelNum,fileDir,ssnDir,varargin)
%% Updated 4/09/2017
%
% This is a script for creating a statMatrix for a single .continuous file
% Example input for running statMatrix script for Channel 1 of a session:
% [statMatrixColIDs,statMatrix, AnalysisColumns] = statMatrixOE(recStart, 100 ,1, fileDir, ssnDir, 'processor', 'ADC', 'CH','session')
% ProcessorNum = 100 (FPGA or Raw); processorNum = 102 (Bandpass)
% fileDir = the folder directory of .continuous, .events, and/or .spike files
% For recStart, open the messages.events file (using Notepad or Notepad++) in target file directory to find starting time of recording
% ssnDir = the folder directory and the file name of your session data

%% Argument String Detection
% Determines what data to extract from data sets

if sum(cell2mat(cellfun(@(a)(strcmp(a,'ADC') | strcmp(a, 'adc')), varargin, 'uniformoutput', 0)))>=1
    [Events_TS,~,~,~] = BehavTS(fileDir,processorNum,recStart,'downsample');
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

%% Bandpass Filter and Hilbert Transform of LFP

% Lowpass Filtered 0-300Hz and Transformed

order = 3;
fcut  = 300;
NyquistFreq = sampleRate/2;
[b,l] = butter(order,fcut/(NyquistFreq),'low');
LFP_Bandpass = filter(b,l,LFP_Raw);

LFPHilb_Bandpass = atan2(imag(hilbert(LFP_Bandpass)), LFP_Bandpass); % Corresponding Hilbert Transform

% Check if signal needs to be downsampled before the Hilbert transforms
if sampleRate > 1000
    n_samples = sampleRate/1000;
    LFP_Bandpass = downsample(LFP_Bandpass,n_samples);
    LFPHilb_Bandpass = downsample(LFPHilb_Bandpass,n_samples);
    TS_Raw = downsample(TS_Raw,n_samples);
end

% Theta filtered and transformed

MinTheta = 4;
MaxTheta = 12;

[LFPHilb_Theta, LFP_Theta] = PhaseFreqDetectAbbr(LFP_Bandpass,TS_Raw,MinTheta,MaxTheta);

% Beta filtered and transformed

MinBeta = 20;
MaxBeta = 40;

[LFPHilb_Beta, LFP_Beta] = PhaseFreqDetectAbbr(LFP_Bandpass,TS_Raw,MinBeta,MaxBeta);

%% Organize LFP Columns
% Organizes all LFP related vectors into columns

TS_Col = transpose(TS_Raw);
LFPBandpass_Col = LFP_Bandpass;
LFPHilbBandpass_Col = LFPHilb_Bandpass;
LFP_Theta_Col = LFP_Theta;
LFPHilb_Theta_Col = LFPHilb_Theta;
LFP_Beta_Col = LFP_Beta;
LFPHilb_Beta_Col = LFPHilb_Beta;


LFP_Mat = [LFPBandpass_Col,LFPHilbBandpass_Col,LFP_Theta_Col,LFPHilb_Theta_Col, LFP_Beta_Col, LFPHilb_Beta_Col];

LFPDataHeaders = ["LFP_Raw","LFP_Raw_Hilbert","LFP_Theta","LFP_Theta_Hilbert","LFP_Beta","LFP_Beta_Hilbert"];


%% Create and Organize Behavioural Columns

% Initializes Behavioral Column Headers for Label References
seqLength = ssnData(1).Settings.SequenceLength;

for seq = 1:seqLength
    behDataHeaders{seq} = ['Odor' num2str(seq)];
    behDataHeaders{seq+seqLength} = ['Position' num2str(seq)];
end

behDataHeaders = [behDataHeaders, "InSeqLog", "Performance"];

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
AllOdorOnset_Index = find(AllOdorOnset);


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
        InSeqLog(AllOdorOnset_Index(SessionLogIndex)) = 0;
    end
end
InSeqLog_Col = transpose(InSeqLog);

% Organizes Position data from Odor Indices
Position1 = histcounts(Odor1_Index,[0 TS_Raw]);
Position2 = histcounts(Odor2_Index,[0 TS_Raw]);
Position3 = histcounts(Odor3_Index,[0 TS_Raw]);
Position4 = histcounts(Odor4_Index,[0 TS_Raw]);

Position1_Col = transpose(Position1);
Position2_Col = transpose(Position2);
Position3_Col = transpose(Position3);
Position4_Col = transpose(Position4);

% Creates Behavrioral Data Matrix 
Behav_Mat = [Odor1_Col,Odor2_Col,Odor3_Col,Odor4_Col,Position1_Col,Position2_Col,...
    Position3_Col,Position4_Col,InSeqLog_Col,Performance_Col];

%% Arrange all matrices/columns into full statMatrix

statMatrix = [TS_Col LFP_Mat Behav_Mat];
statMatrixColIDs = ['TimeBins', LFPDataHeaders, behDataHeaders];


%% Rearrange trial data for analysis
% Extracts indices for In-sequence correct and incorrect, and out-sequence
% correct and incorrect trials for CreateSpectrogram.m

InSeqTrialIndex = find(InSeqLog);
InSeqCorrect = Performance_Col(InSeqTrialIndex) == 1;
InSeqIncorrect = Performance_Col(InSeqTrialIndex) == -1;


InSeqCorrectIndex = InSeqTrialIndex(InSeqCorrect);
InSeqIncorrectIndex = InSeqTrialIndex(InSeqIncorrect);

OutSeqLog = zeros(size(LFP_Bandpass));

for SessionLogIndex = 1:size(SessionLog,2)
    if SessionLog(SessionLogIndex) ~= 0
        OutSeqLog(AllOdorOnset_Index(SessionLogIndex)) = 1;
    else
        OutSeqLog(AllOdorOnset_Index(SessionLogIndex)) = 0;
    end
end

OutSeqTrialIndex = find(OutSeqLog);

OutSeqCorrect = Performance_Col(OutSeqTrialIndex) == 1;
OutSeqIncorrect = Performance_Col(OutSeqTrialIndex) == -1;


OutSeqCorrectIndex = OutSeqTrialIndex(OutSeqCorrect);
OutSeqIncorrectIndex = OutSeqTrialIndex(OutSeqIncorrect);

% I totally forgot why I put these two values to zero
OutSeqCorrectIndex(201) = 0;
OutSeqIncorrectIndex(201) = 0;

AnalysisColumns = [InSeqCorrectIndex',OutSeqCorrectIndex,InSeqIncorrectIndex',OutSeqIncorrectIndex];


end



