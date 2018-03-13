function [Hilbert_phase, signal_filtered] = PhaseFreqDetectAbbr(signal, signal_ts,MinFreq,MaxFreq)

% Outputs
% 1) Hilbert_phase - Phase of oscillation determined by the hilbert transform
% 2) signal_filtered - bandpass filtered signal

% signal = Tchannel#_#(:,1), voltage in vectors; 
% signal_ts%=Tchannel#_#(:,2), timestamp (commented by Damon, Jan 14 2016)

% signal_filtered should be bandpass filtered (filtfilt).  It's most
% important to get high frequency oscillations out - so 1-60 is okay, 1-30
% or even 4-12Hz might be better, depending on your need. 

% MinFreq will set the boundaries of your bandpass filter AND the parameters 
% tolerated for peak and trough detection.  These can be changed to whatever 
% you want.  6-10 is about standard.  This will not filter your signal - it 
% will only constrain this code to find cycles within this range. 

%MinFreq = 20; MaxFreq = 40;

% 10/28/16  Modified Hilbert transform to be performed on the filtered
%           rather than the raw signal - GE

%% Bandpass filter
if nargin > 2
    fs = 1/(signal_ts(2)-signal_ts(1));

    Wn_FRange = [MinFreq/(fs/2) MaxFreq/(fs/2)]; % normalized by the nyquist frequency

    [bFRange, aFRange] = butter(3,Wn_FRange);

    signal_filtered = filtfilt(bFRange,aFRange,signal);
else
    signal_filtered = signal;
end
%% Hilbert transform

Hilbert_phase = atan2(imag(hilbert(signal_filtered)), signal_filtered);
