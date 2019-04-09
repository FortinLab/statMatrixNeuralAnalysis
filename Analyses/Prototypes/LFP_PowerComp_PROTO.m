function LFP_PowerComp(recStart, fileDir, ssnDir, channelNum, processorNum)
%% Determine number of files from target folder

% origDir = cd; % saves current file directory
% cd(fileDir) % changes to target directory
% Filename = sprintf('%d_CH',processorNum);
% numCHfiles = length(dir([Filename,'*']));
% cd(origDir) % reverts to original directory

%% Extract Raw LFP and Z-score
% for channelNum = 1:numCHfiles
windowSize = 4000;

%% Create Stat Matrix
[~,statMatrix,AnalysisColumns] = CreateStatMatrixFromOE(recStart, processorNum, channelNum, fileDir, ssnDir,...
    'processor', 'ADC', 'CH','session');
AnalysisColumns(AnalysisColumns>2472289) = 0;

%% Extract different trial periods
[Events_TS,OnTS,~,Raw_TS] = BehavTS(fileDir,processorNum,recStart,'downsample');
Raw_TS = Raw_TS(1:2472289);
deleteMe = fieldnames(OnTS);
for dM = 1:length(deleteMe)
    curValsToDelete = OnTS.(deleteMe{dM});
    curValsToDelete(curValsToDelete>((2472289/1000)+(windowSize/2/1000))) = [];
    OnTS.(deleteMe{dM}) = curValsToDelete;
end
deleteMeToo = fieldnames(Events_TS);
for dMt = 1:length(deleteMeToo)
    curValsToDeleteToo = Events_TS.(deleteMeToo{dMt});
    Events_TS.(deleteMeToo{dMt}) = curValsToDeleteToo(1:2472289);
end

% Combines all odor trials for later use
AllOdorOnset = Events_TS.OnPeakTS_1 + Events_TS.OnPeakTS_2 +...
    Events_TS.OnPeakTS_3 + Events_TS.OnPeakTS_4;
AllOdorOnset_Index = find(AllOdorOnset==1);
AllOdorOnset_TS = Raw_TS(AllOdorOnset_Index);

% Poke Out Time Periods
for SeqType = 1:2 % 1 = InSeq; 2 = OutSeq
    for C_vs_InC = 1:2 % 1 = Correct; 2 = Incorrect
        if SeqType == 1
            % Poke Out
            if C_vs_InC == 1 % InSeq Correct
                InSeqCorrectIndex = AnalysisColumns(AnalysisColumns(:,1)>0,1);
                InSeqCorrectTS = Raw_TS(InSeqCorrectIndex);
                PokeOutCenters_IC = nearestpoint(InSeqCorrectTS,OnTS.OnPeakTS_6,'next'); % First vector ...
                % is the one positioned before, so the opposite argument must be indicated
                PokeOutCenters_IC_TS = OnTS.OnPeakTS_6(PokeOutCenters_IC);
                [~,~,PokeOutCenters_IC_Index] = intersect(PokeOutCenters_IC_TS,Raw_TS);
            end
            if C_vs_InC == 2 % InSeq Incorrect
                % Poke Out
                InSeqIncorrectIndex = AnalysisColumns(AnalysisColumns(:,3)>0,3);
                InSeqIncorrectTS = Raw_TS(InSeqIncorrectIndex);
                PokeOutCenters_II = nearestpoint(InSeqIncorrectTS,OnTS.OnPeakTS_8,'next');
                PokeOutCenters_II_TS = OnTS.OnPeakTS_8(PokeOutCenters_II);
                [~,~,PokeOutCenters_II_Index] = intersect(PokeOutCenters_II_TS,Raw_TS);
            end
        end
        if SeqType == 2
            if C_vs_InC == 1 % OutSeq Correct
                OutSeqCorrectIndex = AnalysisColumns(AnalysisColumns(:,2)>0,2);
                OutSeqCorrectTS = Raw_TS(OutSeqCorrectIndex);
                PokeOutCenters_OC = nearestpoint(OutSeqCorrectTS,OnTS.OnPeakTS_6,'next');
                PokeOutCenters_OC_TS = OnTS.OnPeakTS_6(PokeOutCenters_OC);
                [~,~,PokeOutCenters_OC_Index] = intersect(PokeOutCenters_OC_TS,Raw_TS);
            end
            if C_vs_InC == 2 % OutSeq Incorrect
                OutSeqIncorrectIndex = AnalysisColumns(AnalysisColumns(:,4)>0,4);
                OutSeqIncorrectTS = Raw_TS(OutSeqIncorrectIndex);
                PokeOutCenters_OI = nearestpoint(OutSeqIncorrectTS,OnTS.OnPeakTS_8,'previous');
                PokeOutCenters_OI_TS = OnTS.OnPeakTS_8(PokeOutCenters_OI);
                [~,~,PokeOutCenters_OI_Index] = intersect(PokeOutCenters_OI_TS,Raw_TS);
            end
        end
    end
end

PokeOutCenters_C_Index = union(PokeOutCenters_IC_Index,PokeOutCenters_OC_Index);
PokeOutCenters_I_Index = union(PokeOutCenters_II_Index,PokeOutCenters_OI_Index);

% 1st Error Gate = Check if there are overlapping time indices for Poke Out
if length(PokeOutCenters_C_Index) == length(PokeOutCenters_IC_Index) + length(PokeOutCenters_OC_Index)
    disp('No overlaps for Poke Out Correct Indices')
else
    disp('ERROR: Overlap detected for Poke Out Correct Indices')
    return
end

if length(PokeOutCenters_I_Index) == length(PokeOutCenters_II_Index) + length(PokeOutCenters_OI_Index)
    disp('No overlaps for Poke Out Incorrect Indices')
else
    disp('ERROR: Overlap detected for Poke Out Incorrect Indices')
    return
end

PokeOutCheck = union(PokeOutCenters_I_Index,PokeOutCenters_C_Index);

if length(PokeOutCheck) == length(PokeOutCenters_I_Index) + length(PokeOutCenters_C_Index)
    disp('No overlaps for Poke Out Indices')
else
    disp('ERROR: Overlap detected for Poke Out Indices')
    return
end

% Poke In Time Periods
% I didn't combine Poke Out timestamps like I did for Poke In in case
% we needed to discern InSec vs. OutSeq. Poke In always happens for the
% 1st odor so that is why they are just combined.

CorrectTS = union(InSeqCorrectTS,OutSeqCorrectTS);
PokeInCenters_C = nearestpoint(CorrectTS,AllOdorOnset_TS,'next');
PokeInCenters_C_TS = AllOdorOnset_TS(PokeInCenters_C);
[~,~,PokeInCenters_C_Index] = intersect(PokeInCenters_C_TS,Raw_TS);

IncorrectTS = union(InSeqIncorrectTS,OutSeqIncorrectTS);
PokeInCenters_I = nearestpoint(IncorrectTS,AllOdorOnset_TS,'next');
PokeInCenters_I_TS = AllOdorOnset_TS(PokeInCenters_I);
[~,~,PokeInCenters_I_Index] = intersect(PokeInCenters_I_TS,Raw_TS);

PokeInCheck = union(PokeInCenters_I_Index,PokeInCenters_C_Index);

% 2nd Error Gate = Check if there are overlapping time indices for
% Poke In
if length(PokeInCheck) == length(PokeInCenters_I_Index) + length(PokeInCenters_C_Index)
    disp('No overlaps for Poke In Indices')
else
    disp('ERROR: Overlap detected for Poke In Indices')
    return
end

% 3rd and Final Error Gate = Check if there are overlaps between Poke in and
% out
PokeFinalCheck = union(PokeInCheck,PokeOutCheck);

if length(PokeFinalCheck) == length(PokeInCheck) + length(PokeOutCheck)
    disp('No overlaps for Poke In and Poke Out Indices')
else
    disp('ERROR: Overlap detected for Poke In and Poke Out Indices')
    return
end
%% Run the analysis


figure
for freqBand = 1:6
    % Frequency Band key
    % Theta = 4-12 Hz
    % Low Beta = 13-19 Hz
    % Beta = 20-40 Hz
    % Low Gamma = 41-59 Hz - Skips 60 Hz because that is the freq for electrical noise
    % High Gamma = 61-80 Hz
    % Ripple = 150-250 Hz
    
    statMatrix = statMatrix(1:2472290,:);
    if freqBand == 1
        MinFreq = 4;
        MaxFreq = 12;
        wSize = 1.25*(MaxFreq+MinFreq)/2;
        plotLabel = 'Theta';
    end
    if freqBand == 2
        MinFreq = 13;
        MaxFreq = 19;
        wSize = 1.25*(MaxFreq+MinFreq)/2;
        plotLabel = 'Low Beta';
    end
    if freqBand == 3
        MinFreq = 20;
        MaxFreq = 40;
        wSize = 1.25*(MaxFreq+MinFreq)/2;
        plotLabel = 'Beta';
    end
    if freqBand == 4
        MinFreq = 41;
        MaxFreq = 59;
        wSize = 1.25*(MaxFreq+MinFreq)/2;
        plotLabel = 'Low Gamma';
    end
    if freqBand == 5
        MinFreq = 61;
        MaxFreq = 80;
        wSize = 1.25*(MaxFreq+MinFreq)/2;
        plotLabel = 'High Gamma';
    end
    if freqBand == 6
        MinFreq = 150;
        MaxFreq = 250;
        wSize = 1.25*(MaxFreq+MinFreq)/2;
        plotLabel = 'Ripple';
    end
    wSize = floor(wSize);
    [~, LFP_BandFiltered] = PhaseFreqDetectAbbr(statMatrix(:,2),statMatrix(:,1),MinFreq,MaxFreq);
    
    %     windowNdxs = (1/wSize)*1000;
    %     hlfWndo = windowNdxs/2;
    %     paddedVolt = [nan(hlfWndo,1); LFP_BandFiltered; nan(hlfWndo,1)];
    %     tempCurRMS = nan(size(LFP_BandFiltered));
    %     parfor ndx = 1:length(LFP_BandFiltered)
    %         tempCurRMS(ndx) = rms(paddedVolt(ndx:ndx+windowNdxs), 'omitnan');
    %     end
    %     LFP_rms = tempCurRMS;
    
    [LFP_rms,~] = mvrms_mfile_simulink(statMatrix(:,1),LFP_BandFiltered,wSize,0);
    
    
    LFP_BandFiltered_zNorm = zscore(LFP_rms);
    
    %%
    
    
    PokeIn_C_TPs = nan(length(PokeInCenters_C_Index), windowSize+1);
    for numTPs = 1:length(PokeInCenters_C_Index)
        PokeIn_C_TPs(numTPs,:) = LFP_BandFiltered_zNorm(PokeInCenters_C_Index(numTPs)...
            -(windowSize/2):PokeInCenters_C_Index(numTPs)+(windowSize/2));
    end
    
    PokeIn_I_TPs = nan(length(PokeInCenters_I_Index), windowSize+1);
    for numTPs = 1:length(PokeInCenters_I_Index)
        PokeIn_I_TPs(numTPs,:) = LFP_BandFiltered_zNorm(PokeInCenters_I_Index(numTPs)...
            -(windowSize/2):PokeInCenters_I_Index(numTPs)+(windowSize/2));
    end
    
    PokeOut_C_TPs = nan(length(PokeOutCenters_C_Index), windowSize+1);
    for numTPs = 1:length(PokeOutCenters_C_Index)
        PokeOut_C_TPs(numTPs,:) = LFP_BandFiltered_zNorm(PokeOutCenters_C_Index(numTPs)...
            -(windowSize/2):PokeOutCenters_C_Index(numTPs)+(windowSize/2));
    end
    
    PokeOut_I_TPs = nan(length(PokeOutCenters_I_Index), windowSize+1);
    for numTPs = 1:length(PokeOutCenters_I_Index)
        PokeOut_I_TPs(numTPs,:) = LFP_BandFiltered_zNorm(PokeOutCenters_I_Index(numTPs)...
            -(windowSize/2):PokeOutCenters_I_Index(numTPs)+(windowSize/2));
    end
    
    PokeIn_C_meanTP = mean(PokeIn_C_TPs);
    PokeIn_I_meanTP = mean(PokeIn_I_TPs);
    PokeOut_C_meanTP = mean(PokeOut_C_TPs);
    PokeOut_I_meanTP = mean(PokeOut_I_TPs);
    
    time = (-1*windowSize/2/1000):0.001:(windowSize/2/1000);
    
    piC = subplot(2,2,1);
    plot(time,PokeIn_C_meanTP)
    hold on
    title('Poke In - Correct')
    xlabel('Time (sec)')
    ylabel('Z-Normalized Power')
    
    
    poC = subplot(2,2,2);
    plot(time,PokeOut_C_meanTP)
    hold on
    title('Poke Out - Correct')
    xlabel('Time (sec)')
    ylabel('Z-Normalized Power')
    
    piI = subplot(2,2,3);
    plot(time,PokeIn_I_meanTP)
    hold on
    title('Poke In - Incorrect')
    xlabel('Time (sec)')
    ylabel('Z-Normalized Power')
    
    poI = subplot(2,2,4);
    plot(time,PokeOut_I_meanTP)
    hold on
    title('Poke Out - Incorrect')
    xlabel('Time (sec)')
    ylabel('Z-Normalized Power')
    drawnow;
    
    legendLabels{freqBand} = plotLabel;
end
legend(legendLabels{1},legendLabels{2},legendLabels{3},legendLabels{4},legendLabels{5},...
    legendLabels{6})
linkaxes([piC, poC, piI, poI], 'xy');

% folder = strcat('C:\Users\Gabo\Documents\MATLAB\Glenn_spec_Aged 1-17_Session204_ZNormPower\'); %Insert target folder here
%
% if processorNum == 100
%     pngFileName = sprintf('Channel%d_Tetrode%d_FPGA.png', channelNum,tet); % Use for spectrograms in the FPGA node; processor = 100
% elseif processorNum == 102
%     pngFileName = sprintf('Channel%d_Tetrode%d_BP.png', channelNum,tet); % Use for spectrograms in the bandpass node; processor = 102
% end
% fullFileName = fullfile(folder, pngFileName);
%
% saveas(gcf,fullFileName)
% toc
% end
end