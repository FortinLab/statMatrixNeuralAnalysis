
origDir = cd;
[fileDir] = uigetdir(origDir);
if fileDir==0
    disp('Analysis Cancelled')
    return
else
    cd(fileDir)
end

%% Define Parameters
alignment = 'PokeIn';
windowStart = -0.8;
windowEnd = 1.5;
dataBinSize = 125;
freqStart = 2;
freqEnd = 100;
%% Load Relevant Data
% 'trialInfo' : Matrix containing information about the identity and
%           outcome of each trial. Rows represent individual trials;
%           columns are organized thusly:
%               Column 1) Trial performance. 1 = Correct, 0 = Incorrect
%               Column 2) Trial InSeq log. 1 = InSeq, 0 = OutSeq
%               Column 3) Trial Position.
%               Column 4) Trial Odor.
[~, ~, lfpEpoch, lfpIDs, ~, eventTimeBins, trialInfo] = EpochExtraction_SM(alignment, windowStart, windowEnd, 'org', 'TiUTr');

%% 
lfpVals = repmat({nan(length(freqStart:freqEnd),length(windowStart:0.001:windowEnd))}, [1,length(lfpIDs),size(lfpEpoch,3)]);
% lfpVals = cell(1,length(lfpIDs), size(lfpEpoch,3));
for tet = 1:length(lfpIDs)
    fprintf('Starting tetrode %s...', lfpIDs{tet});
    for trl = 1:size(lfpEpoch,3)
        lfpVals{1,tet,trl} = MorletAG(lfpEpoch(:,tet,trl), round(1/mode(diff(eventTimeBins))), freqStart, freqEnd);
    end
    figure;
    imagesc(windowStart:0.001:windowEnd, freqStart:freqEnd,mean(cell2mat(lfpVals(1,tet,:)),3))
    set(gca, 'ydir', 'normal')
    line(get(gca, 'xlim'), [3.5 3.5], 'color', 'white', 'linewidth', 2, 'linestyle', '--')
    line(get(gca, 'xlim'), [11.5 11.5], 'color', 'white', 'linewidth', 2, 'linestyle', ':')
    line(get(gca, 'xlim'), [7.5 7.5], 'color', 'white', 'linewidth', 2, 'linestyle', '-.')
    line(get(gca, 'xlim'), [31.5 31.5], 'color', 'white', 'linewidth', 2)
    line(get(gca, 'xlim'), [15.5 15.5], 'color', 'white', 'linewidth', 2)
    title(lfpIDs{tet}, 'interpreter', 'none');
    drawnow
    fprintf('Completed\n');
end

    
