function CreateSpectrogram(processorNum,channelNum,recStart,ssnDir,fileDir)
%% Updated 3/07/2017
%% Update Changes: Changed script so that time window for Morlet is inputted by the user
% This is a script for creating spectrograms for each channel
% Example input for running CreateSpectrogram script for Channel 1 of a session:
% CreateSpectrogram(processorNum,1,recStart,ssnDir,fileDir,[minTimestamp:1/SamplingFrequency:maxTimestamp);
% NOTE: ssnDir includes the file name as well, fileDir is only the folder
% directory of the files
% Let Sampling Frequency = 1000
% Normal Setting: time = -1.5s : 1/SamplingFrequency: 1.5s
% For a rundown on what the inputs are, read comments on statMatrixOE.m

%% Extract statMatrix Data
[~,statMatrix,AnalysisColumns] = CreateStatMatrixFromOE(recStart, processorNum, channelNum, fileDir, ssnDir,...
    'processor', 'ADC', 'CH','session');

InSeqCorrect = AnalysisColumns(:,1);
OutSeqCorrect = AnalysisColumns(:,2);

InSeqCorrectIndex = InSeqCorrect(InSeqCorrect>0);
OutSeqCorrectIndex = OutSeqCorrect(OutSeqCorrect>0);

lfp = statMatrix(:,2);

%% Wavelet Parameters       
Fs = 1000;                                                  %% Sampling Rate in Hertz
time = -1.5:1/Fs:1.5;                                       %% Length of Time Wavelet extends over
min_freq =  2;                                              %% Smallest Morlet
max_freq = 100;                                             %% Largest Morlet
freq = min_freq:max_freq;                                   %% Frequency of Morlet Wavelet
% time = [-1.5,1/1000,1.5];



%% Butterworth Notch Filter at 60hz
[b1, a1] = butter(2, [.118 .122], 'stop');                   %% Remove n+1 hz (noise); [0.118 0.122] is the (fc-1Hz)/(fs/2) and (fc+1Hz)/(fs/2)
avgM = zeros(length(freq), length(time), 2);                %% Store the Averages for InSeq and OutSeq for wire
tic
for SequenceType = 1:2
        if SequenceType == 1
            LFPIndex = InSeqCorrectIndex;
            TrialMatrix = zeros(length(freq), length(time), []);  
        elseif SequenceType == 2
            LFPIndex = OutSeqCorrectIndex;
            TrialMatrix = zeros(length(freq), length(time), []);  
        end     
    for i = 1:length(LFPIndex)

        
        data = lfp(LFPIndex(i)-1500: LFPIndex(i)+1500);   %% Epoch - 1.5 seconds in both directions (3secs)
        data = data';                                           %% Transpose to 1 row of columns.
        data_filt = filtfilt(b1, a1, data);                     %% Apply Filter
        
        %% Create Complex Wave w/Variable Frequency
        for fi = 1:length(freq)
            sin_wave = exp(1i*2*pi*freq(fi).*time);             %% Complex Sine Wave
            
            s = 7/(2*pi*freq(fi));                              %% Standard Dev/ Num Cycles for Gaussian
            gaus = exp((-time.^2) ./ (2*s^2));                  %% Complex Gaussian
            
            morlet = sin_wave .* gaus;                          %% multiply Sine and Gaussian
            
            len_data = length(data_filt);                       %% Length of Data
            len_kernal = length(morlet);                        %% Length of Kernal
            len_conv = (len_data + len_kernal - 1);             %% Length of Convoluted Data
            half = floor(len_kernal/2)+1;                       %% Used when matching analytic signal length
            
            dataX = fft(data_filt, len_conv);                   %% Fast Fourier Transform of Data (correct length)
            
            morletX = fft(morlet,len_conv);                     %% Fast Fourier Transform of Morlet
            morletX = morletX ./ max(morletX);                  %% Normalize to largest
            
            conv_data =  (dataX .* morletX);                    %% Convolve - Element-wise mutliplication
            
            analytic_signal = ifft(conv_data);                  %% Return to Time-Domain
            analytic_signal = analytic_signal(half-1:end-half); %% Match Length
            
            power = abs(analytic_signal).^2;                    %% Extract Power
            power = (10*log10(power));                          %% Convert to Logarithmic Scale
            TrialMatrix(fi, : , i) = power;                     %% Add Data to 3-Dimensional Matrix
        end
    end
    %% Average Trial
    if SequenceType == 1
        avgM( : , : , 1) = mean(TrialMatrix, 3);                      %% Average All Trials for Current Tetrode (InSeq)
    elseif SequenceType == 2
        avgM( : , : , 2) = mean(TrialMatrix, 3);                      %% Average All Trials for Current Tetrode (OutSeq)
    end
end
toc
%% Time-Frequency Plot

if channelNum <= 4
    tet = 9;
elseif channelNum > 4 && channelNum <=8
    tet = 10;
elseif channelNum > 8 && channelNum <=12
    tet = 11;
elseif channelNum > 12 && channelNum <=16
    tet = 12;
elseif channelNum > 16 && channelNum <=20
    tet = 0;
elseif channelNum > 20 && channelNum <=24
    tet = 0;
elseif channelNum > 24 && channelNum <=28
    tet = 15;
elseif channelNum > 28 && channelNum <=32
    tet = 16;
elseif channelNum > 32 && channelNum <=36
    tet = 1;
elseif channelNum > 36 && channelNum <=40
    tet = 2;
elseif channelNum > 40 && channelNum <=44
    tet = 3;
elseif channelNum > 44 && channelNum <=48
    tet = 4;
elseif channelNum > 48 && channelNum <=52
    tet = 0;
elseif channelNum > 52 && channelNum <=56
    tet = 0;
elseif channelNum > 56 && channelNum <=60
    tet = 0;
elseif channelNum > 60 && channelNum <=64
    tet = 8;
end

for jj = 1:2                                                %% In Seq vs. Out Seq
    if jj == 1                                              %% Info for Subplot
        type = 'In Sequence';
        SubPlotPosition = 1;
        SequenceType = 1;
    else
        type = 'Out Sequence';
        SequenceType = SequenceType + 1;
        SubPlotPosition = SubPlotPosition +1;
    end
    
    figure(1)
    subplot(2,1,SubPlotPosition)
    contourf(time, freq, avgM( : , : , SequenceType),...        %% Create 3D Graph
        40, 'linecolor', 'none')
    title(sprintf('Tetrode %d - %s',...                         %% Graph Title
        tet, type), 'FontSize', 14)                  %% Graph Font
    xlabel('Time (seconds)', 'FontSize', 16),...                %% Graph X Label
        ylabel('Frequency (hz)', 'FontSize', 16)                    %% Graph Y Label
    %     caxis([-38 -18])                                          %% Color Ax
    colormap(jet)                                               %% Color scheme for Colorbar
    colorbar                                                    %% Open Colorbar
end
end

