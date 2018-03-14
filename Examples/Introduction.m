%% Using the StatMatrix Example %%

% The easiest way to extract information from the statMatrix is to use 
% logical operations. Below are a few examples of how you might use this to
% work with the data in the statMatrix Format.

% First, let's load a data file from SuperChris' dataset. Make sure you're
% in the directory containing SuperChris' data, or you've added that 
% directory to your path.
load('SuperChris_02-12-09_MS_T14.mat');

%% Extracting Unit Timestamps
% The first thing to do is identify the column associated with your unit of
% interest. In this situation you know you want to look at unit T14-U4.

% First you'll need to identify it's location within the statMatrix.

t14u4ColPosLog = strcmp(statMatrixColIDs, 'T14-U4');

% Now, with your logical row in hand you can extract the spiking
% activity from the statMatrix and create a logical column indicating when
% the spikes occurred.

t14u4SpikeLog = statMatrix(:,t14u4ColPosLog)==1;

% We can now use the logical column to extract the timestamps associated
% with each spike. Remember, the first column in the statMatrix is always
% the timestamp vector.

t14u4TSs = statMatrix(t14u4SpikeLog, 1);

%% Extracting Trial Timestamps
% In this case, let's extract the timestamps for every time odor 2 (B) is
% presented and the animal got it correct and separate out InSeq and
% OutSeq presentations.

% First, create your logical row to identify the odor 2 column.
odor2ColLog = strcmp(statMatrixColIDs, 'Odor2');
% Next, create a logical column to identify when odor 2 was presented.
odor2PresLog = statMatrix(:,odor2ColLog);

% Now identify the performance column and create a performance logical
% column
perfLogColLog = strcmp(statMatrixColIDs, 'PerformanceLog');
perfLog = statMatrix(:,perfLogColLog)==1;
% ... and use it to select only odor 2 trials he got correct.
odor2PresLog = odor2PresLog & perfLog;

% Now that we have a logical column of odor 2 presentations, let's create a
% logical column to identify the InSeq trials.
inSeqLogColLog = strcmp(statMatrixColIDs, 'InSeqLog');
inSeqTrlLog = statMatrix(:,inSeqLogColLog)==1;
outSeqTrlLog = statMatrix(:,inSeqLogColLog)==-1;

% Now, with these three logicals, odor2PresLog, inSeqTrlLog, outSeqTrlLog
% we have everything we would need to determine when odor 2 was presented
% InSeq and when it was presented OutSeq.
odor2InSeqTSs = statMatrix(odor2PresLog & inSeqTrlLog, 1);
odor2OutSeqTSs = statMatrix(odor2PresLog & outSeqTrlLog, 1);

%% Extract Unit Trial Activity
% Now let's combine the two. With out unit and trial timestamps we can now
% examine how the units were active during the different trial periods.

% To do this we need to use a for loop to extract activity during a period
% around trial start. Let's first define our period of interest. Let's
% look at 500ms before the trial starts and the first half of the trial,
% 500ms. Remember, the timestamps are in seconds.
trialWindow = [-0.5 0.5];

% Now let's go through each trial and extract the spiking activity relative
% to that trial start. Let's do the InSeq trials first.
% First create an empty cell vector, this is where we're going to put the
% data we extract.
inSeqTrlSpkTS = cell(length(odor2InSeqTSs),1);
for trl = 1:length(odor2InSeqTSs)
    % First identify the current trial start time
    curTrlStart = odor2InSeqTSs(trl);
    % Now let's align all our spike timestamps to trial start by
    % subtracting the two.
    spikesRelativeToCurTrlStart = t14u4TSs - curTrlStart;
    % Now create a logical vector to select only those spikes that occur in
    % our window of interest.
    curTrlWindowSpksLog = spikesRelativeToCurTrlStart>trialWindow(1) & spikesRelativeToCurTrlStart<trialWindow(2);
    % Now use this logical vector to fill in our cell array data we created
    % outside
    inSeqTrlSpkTS{trl} = spikesRelativeToCurTrlStart(curTrlWindowSpksLog);
end

% Now let's do the same thing with the OutSeq Trials
outSeqTrlSpkTS = cell(length(odor2OutSeqTSs),1);
for trl = 1:length(odor2OutSeqTSs)
    curTrlStart = odor2OutSeqTSs(trl);    
    spikesRelativeToCurTrlStart = t14u4TSs - curTrlStart;
    curTrlWindowSpksLog = spikesRelativeToCurTrlStart>trialWindow(1) & spikesRelativeToCurTrlStart<trialWindow(2);
    outSeqTrlSpkTS{trl} = spikesRelativeToCurTrlStart(curTrlWindowSpksLog);
end

% With the data organized this way it's easy to now combine them all and
% make a simple PSTH and make separate plots for InSeq and OutSeq 
figure;
inSeqPSTH = subplot(1,2,1);
% Use cell2mat to take all your trial oriented spike times and then use
% histcounts to give you tallies of how many fall into each time bin
inSeqCounts = histcounts(cell2mat(inSeqTrlSpkTS), trialWindow(1):0.05:trialWindow(2));
% Divide by the number of trials to enable comparisons between InSeq and
% OutSeq conditions
inSeqCounts = inSeqCounts/length(odor2InSeqTSs);
histogram('BinEdges', trialWindow(1):0.05:trialWindow(2), 'BinCounts', inSeqCounts);
outSeqPSTH = subplot(1,2,2);
outSeqCounts = histcounts(cell2mat(outSeqTrlSpkTS), trialWindow(1,:):0.05:trialWindow(2));
outSeqCounts = outSeqCounts/length(odor2OutSeqTSs);
histogram('BinEdges', trialWindow(1):0.05:trialWindow(2), 'BinCounts', outSeqCounts);
linkaxes([inSeqPSTH, outSeqPSTH], 'xy')

