% function CreateSpectrogram_Aaron(channelNumber,recStart,fileDir)


% [~,statMatrix,AnalysisColumns] = statMatrixOE(recStart, channelNumber , fileDir, 'processor', 'ADC', 'CH','session');
%
% InSeqCorrect = AnalysisColumns(:,1);
% OutSeqCorrect = AnalysisColumns(:,2);
%
% InSeqCorrectIndex = InSeqCorrect(InSeqCorrect>0);
% OutSeqCorrectIndex = OutSeqCorrect(OutSeqCorrect>0);
%
% lfp = statMatrix(:,2);



%% Wavelet Parameters
Fs = 1000;                                                  %% Sampling Rate in Hertz
time = -1.5:1/Fs:1.5;                                       %% Length of Time Wavelet extends over
min_freq =  2;                                              %% Smallest Morlet
max_freq = 100;                                             %% Largest Morlet
freq = min_freq:max_freq;                                   %% Frequency of Morlet Wavelet



%% Butterworth Notch Filter at 60hz
[b1, a1] = butter(2, [.118 .122], 'stop');                   %% Remove n+1 hz (noise)
avgM = zeros(length(freq), length(time), 24);                %% Store the Averages for InSeq and OutSeq for wire

%% Indexes
InSeqCorrectIndex = [72, 73, 78, 82, 94, 95, 101, 102, 103, 106, 107, 116, 126,...
    133, 138, 144, 160, 162, 168, 174, 178, 179, 182, 190];
InSeqCorrectIndex = InSeqCorrectIndex';

tetrodes = [2,9, 8, 7, 15, 13, 14, 23, 16, 22, 19, 20];

for ci = 1:2
    
    for ti = 1:length(tetrodes)
        TrialMatrix = zeros(length(freq), length(time), 25);
        
        for ri = 1:length(InSeqCorrectIndex)
            if ci == 1                                              %% Choose Behavior based on Column
                row_number = indices(InSeqCorrectIndex(ri), ci);                %% Find Row Number Reference from Indices Variable
            elseif ci == 2
                row_number = indices(ri, ci);                       %% Find Row Number Reference from Indices Variable
            elseif ci == 3
                row_number = sd_pos(ri, ci);                        %% Find Row Number Reference from Position Variable
            end
            
            a = behavior(row_number, 6);                            %% Time 0 to be taken from Sample
            data = lfp(a-1500: a+1500, ti);                         %% Epoch - 1.5 seconds in both directions (3secs)
            data = data';                                           %% Transpose to 1 row of columns.
            data_filt = filtfilt(b1, a1, data);
            
            
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
                TrialMatrix(fi, : , ri) = power;                     %% Add Data to 3-Dimensional Matrix
            end
            
            
        end
        %% Average Trial
        if ci == 1
            avgM( : , : , ti) = mean(TrialMatrix, 3);                      %% Average All Trials for Current Tetrode (InSeq)
        elseif ci == 2
            avgM( : , : , ti+12) = mean(TrialMatrix, 3);                      %% Average All Trials for Current Tetrode (OutSeq)
        end
    end
end


%% Time-Frequency Plot
tic

% if channelNumber <= 4
%     tet = 14;
% elseif channelNumber > 4 && channelNumber <=8
%     tet = 13;
% elseif channelNumber > 8 && channelNumber <=12
%     tet = 9;
% elseif channelNumber > 12 && channelNumber <=16
%     tet = 10;
% elseif channelNumber > 16 && channelNumber <=20
%     tet = 11;
% elseif channelNumber > 20 && channelNumber <=24
%     tet = 12;
% elseif channelNumber > 24 && channelNumber <=28
%     tet = 16;
% elseif channelNumber > 28 && channelNumber <=32
%     tet = 15;
% elseif channelNumber > 32 && channelNumber <=36
%     tet = 22;
% elseif channelNumber > 36 && channelNumber <=40
%     tet = 21;
% elseif channelNumber > 40 && channelNumber <=44
%     tet = 17;
% elseif channelNumber > 44 && channelNumber <=48
%     tet = 18;
% elseif channelNumber > 48 && channelNumber <=52
%     tet = 19;
% elseif channelNumber > 52 && channelNumber <=56
%     tet = 20;
% elseif channelNumber > 56 && channelNumber <=60
%     tet = 24;
% elseif channelNumber > 60 && channelNumber <=64
%     tet = 23;
% end
a = 1;                                                              %% Data from Matrix Sheet( , , a)
b = 1;                                                              %% Position within SubPlots
tet = 1;                                                            %% Tetrode Number
for kk = 1:4                                                        %% Create 4 Figures
    for jj = 1:2                                                    %% In Seq vs. Out Seq
        for ii = 1:3                                                %% SubPlot Column Number
            if jj == 1                                              %% Info for Subplot
                type = 'In Sequence';
                b = ii;
                if kk == 1
                    a = ii;
                    tet = ii;
                elseif kk == 2
                    a = (ii + 3);
                    tet = (ii + 3);
                elseif kk == 3
                    a = (ii + 6);
                    tet = (ii + 6);
                else
                    a = (ii + 9);
                    tet = (ii + 9);
                end
            else
                type = 'Out Sequence';
                b = (ii + 3);
                if kk == 1
                    a = (ii + 12);
                    tet = (ii);
                elseif kk == 2
                    a = (ii + 15);
                    tet = (ii + 3);
                elseif kk == 3
                    a = (ii + 18);
                    tet = (ii + 6);
                else
                    a = (ii + 21);
                    tet = (ii + 9);
                end
            end
            
            figure(kk)
            subplot(2, 3, b)
            contourf(time, freq, avgM( : , : , a),...                   %% Create 3D Graph
                40, 'linecolor', 'none')
            title(sprintf('Tetrode_%d - %s',...                         %% Graph Title
                tetrodes(tet), type), 'FontSize', 14)                  %% Graph Font
            xlabel('Time (seconds)', 'FontSize', 16),...                %% Graph X Label
                ylabel('Frequency (hz)', 'FontSize', 16)                    %% Graph Y Label
            %     caxis([-38 -18])                                          %% Color Ax
            colormap(jet)                                               %% Color scheme for Colorbar
            colorbar                                                    %% Open Colorbar
        end
    end
end
%end
