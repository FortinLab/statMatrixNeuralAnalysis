function [freqPower] = MorletAG(lfp, sampleRate, minFreq, maxFreq)
%% MorletAG
%   Functionalized version of Analysis_Morlet_Wavelet.m created by Aaron
%   Gudmundson. 
%
%   Inputs:
%       - lfp: Local field potential voltage values over time. Data should
%           be a vector. Use the RAW trace for this, not any of the
%           bandpass filtered ones
%       - sampleRate: Rate at which LFP recordings were taken OR
%           downsampled to. Value is in Hz.
%       
%% Define empty variables if any
if nargin==0
    disp('You dummy you didn''t give in any data!!!');
    return
elseif nargin==1
    disp('You dummy you didn''t give me the samplerate!!!');
elseif nargin==2
    minFreq = 2;
    maxFreq = 120;
end
%% Define wavelet parameters
time = ((1:length(lfp))-ceil(length(lfp)/2))/sampleRate;                                %% Length of Time Wavelet extends over
freq = minFreq:maxFreq;                                           %% Frequency of Morlet Wavelet

%% Bandstop 60Hz signal
[b1, a1] = butter(2, [59/500 61/500], 'stop');                      %% Remove 60hz Harmonic (noise)
bsData = filtfilt(b1, a1, lfp);                                     %% Apply Filter
% bsData = lfp;
freqPower = nan(length(freq),length(time));
for fi = 1:length(freq)                                 %% Create Complex Wave w/Variable Frequency
    sin_wave = exp(1i*2*pi*freq(fi).*time);             %% Complex Sine Wave
    
    s = 7/(2*pi*freq(fi));                              %% Standard Dev/ Num Cycles for Gaussian
    gaus = exp((-time.^2) ./ (2*s^2));                  %% Complex Gaussian
    
    morlet = sin_wave .* gaus;                          %% multiply Sine and Gaussian
    
    len_data = length(bsData);                          %% Length of Data
    len_kernal = length(morlet);                        %% Length of Kernal
    len_conv = (len_data + len_kernal - 1);             %% Length of Convoluted Data
    half = floor(len_kernal/2)+1;                       %% Used when matching analytic signal length
    
    dataX = fft(bsData', len_conv);                      %% Fast Fourier Transform of Data (correct length)
    
    morletX = fft(morlet,len_conv);                     %% Fast Fourier Transform of Morlet
    morletX = morletX ./ max(morletX);                  %% Normalize to largest
    
    conv_data =  (dataX .* morletX);                    %% Convolve - Element-wise mutliplication
    
    analytic_signal = ifft(conv_data);                  %% Return to Time-Domain
    analytic_signal = analytic_signal(half-1:end-half); %% Match Length
    
    power = abs(analytic_signal).^2;                    %% Extract Power
    power = (10*log10(power));                          %% Convert to Logarithmic Scale
    freqPower(fi, :) = power;                                   %% Add Data to 3-Dimensional Matrix
 
%     fprintf('%.2f Prcnt Completed\n', freq/length(moreFreqs));
end
