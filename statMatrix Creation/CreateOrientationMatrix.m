function CreateOrientationMatrix

origDir = cd;
[fileDir] = uigetdir(origDir);
if fileDir==0
    disp('Analysis Cancelled')
    return
else
    cd(fileDir)
end
dirContents = dir(fileDir);
fileNames = {dirContents.name};

%% Load BehaviorMatrix and Orientation_Data Files
load(fileNames{cellfun(@(a)~isempty(a), strfind(fileNames, 'BehaviorMatrix'))});
load(fileNames{cellfun(@(a)~isempty(a), strfind(fileNames, 'Orientation_Data'))});

%% Compare the timestamps from the behavMatrix and orientation_data files.
%   Since the timestamps used in the statMatrix files is constructed based
%   on the LFP sample rate they don't match up identically with the frame
%   indices taken from the .AVI file. Therefore we are expecting some
%   discrepancy to exist between the video timestamps ('FrameTimestamp'
%   column) and the behavMatrix timestamps ('TimeBin' column), but that 
%   discrepancy should be very small.
frameIndices = [orientData.FileIndices{:,strcmp(orientData.FileIndicesColIDs, 'FrameIndex')}];
frameTimestamps = [orientData.FileIndices{:,strcmp(orientData.FileIndicesColIDs, 'FrameTimestamp')}];

% To assess whether they are from the same data we first check to make sure
% the indices expected in the orientData file can be found within the
% behavMatrix. Then we compare the mean and median for the difference
% between the orientData and behavMatrix and as long as it's below a
% floating point threshold odds are they're from the same data and we can
% assume the frameIndices are correct.
if max(frameIndices) > size(behavMatrix,1)
    error('More frameIndices than actual indices in the behavMatrix. Double check files');
elseif mean(behavMatrix(frameIndices,1) - frameTimestamps') - median(behavMatrix(frameIndices,1) - frameTimestamps') >= 1.0e-10
    error('Difference between mean/median for behavior and orientation timestamp differences is larget than 1.0e-10. Double check files or change threhsold for match violations');
end

%% Create the orientMatrix and orientMatrixColIDs variables
orientMatrixColIDs = [{'TimeBin'}, {'PortX'}, {'PortY'}, {'PortAngle'},...
    {'HeadX'}, {'HeadY'}, {'HeadAngle'},...
    {'TailX'}, {'TailY'}, {'TailAngle'},...
    {'HeadTailLength'}, {'HeadPortLength'}, {'PortTailLength'}];
orientMatrix = [behavMatrix(:,1), nan(size(behavMatrix,1), 12)];

%% Fill in the orientation values from the Orientation_Data file
orientMatrix(frameIndices,strcmp(orientMatrixColIDs, 'PortX')) = [orientData.FileIndices{:,strcmp(orientData.FileIndicesColIDs, 'PortX')}];
orientMatrix(frameIndices,strcmp(orientMatrixColIDs, 'PortY')) = [orientData.FileIndices{:,strcmp(orientData.FileIndicesColIDs, 'PortY')}];
orientMatrix(frameIndices,strcmp(orientMatrixColIDs, 'HeadX')) = [orientData.FileIndices{:,strcmp(orientData.FileIndicesColIDs, 'HeadX')}];
orientMatrix(frameIndices,strcmp(orientMatrixColIDs, 'HeadY')) = [orientData.FileIndices{:,strcmp(orientData.FileIndicesColIDs, 'HeadY')}];
orientMatrix(frameIndices,strcmp(orientMatrixColIDs, 'TailX')) = [orientData.FileIndices{:,strcmp(orientData.FileIndicesColIDs, 'TailX')}];
orientMatrix(frameIndices,strcmp(orientMatrixColIDs, 'TailY')) = [orientData.FileIndices{:,strcmp(orientData.FileIndicesColIDs, 'TailY')}];

%% Calculate the angle and length values
posIndices = find(~isnan(orientMatrix(:,2)));

for pos = 1:length(posIndices)
    curPosNdx = posIndices(pos);
    curPortX = orientMatrix(curPosNdx, strcmp(orientMatrixColIDs, 'PortX'));
    curPortY = orientMatrix(curPosNdx, strcmp(orientMatrixColIDs, 'PortY'));
    curHeadX = orientMatrix(curPosNdx, strcmp(orientMatrixColIDs, 'HeadX'));
    curHeadY = orientMatrix(curPosNdx, strcmp(orientMatrixColIDs, 'HeadY'));
    curTailX = orientMatrix(curPosNdx, strcmp(orientMatrixColIDs, 'TailX'));
    curTailY = orientMatrix(curPosNdx, strcmp(orientMatrixColIDs, 'TailY'));
    
    htVal = sqrt((curTailX - curHeadX)^2 + (curTailY - curHeadY)^2);
    orientMatrix(curPosNdx, strcmp(orientMatrixColIDs, 'HeadTailLength')) = htVal;
    hpVal = sqrt((curPortX - curHeadX)^2 + (curPortY - curHeadY)^2);
    orientMatrix(curPosNdx, strcmp(orientMatrixColIDs, 'HeadPortLength')) = hpVal;
    ptVal = sqrt((curTailX - curPortX)^2 + (curTailY - curPortY)^2);
    orientMatrix(curPosNdx, strcmp(orientMatrixColIDs, 'PortTailLength')) = ptVal;
    
    orientMatrix(curPosNdx, strcmp(orientMatrixColIDs, 'PortAngle')) = rad2deg(acos((ptVal^2 + hpVal^2 - htVal^2)/(2*ptVal*hpVal)));
    orientMatrix(curPosNdx, strcmp(orientMatrixColIDs, 'HeadAngle')) = rad2deg(acos((htVal^2 + hpVal^2 - ptVal^2)/(2*htVal*hpVal)));
    orientMatrix(curPosNdx, strcmp(orientMatrixColIDs, 'TailAngle')) = rad2deg(acos((htVal^2 + ptVal^2 - hpVal^2)/(2*htVal*ptVal)));    
end

%% Save the orientationMatrix File
fileParts = strsplit(fileNames{cellfun(@(a)~isempty(a), strfind(fileNames, 'BehaviorMatrix'))}, '_');
save(sprintf('%s_%s_OrientationMatrix', fileParts{1}, fileParts{2}), 'orientMatrix', 'orientMatrixColIDs');

%% Create a summary figure
figure;
scatAll = subplot(2,2,1);
scatter(orientMatrix(:,strcmp(orientMatrixColIDs, 'HeadX')), orientMatrix(:,strcmp(orientMatrixColIDs, 'HeadY')), 15, 'markeredgecolor', 'none', 'markerfacecolor', 'b', 'markerfacealpha', 0.25);
hold on;
scatter(orientMatrix(:,strcmp(orientMatrixColIDs, 'TailX')), orientMatrix(:,strcmp(orientMatrixColIDs, 'TailY')), 15, 'markeredgecolor', 'none', 'markerfacecolor', 'r', 'markerfacealpha', 0.25);
scatter(orientMatrix(:,strcmp(orientMatrixColIDs, 'PortX')), orientMatrix(:,strcmp(orientMatrixColIDs, 'PortY')), 100, 'markeredgecolor', 'k', 'markerfacecolor', 'k')

scatter(nanmean(orientMatrix(:,strcmp(orientMatrixColIDs, 'HeadX'))), nanmean(orientMatrix(:,strcmp(orientMatrixColIDs, 'HeadY'))), 100, 'x', 'markeredgecolor', 'k', 'markerfacecolor', 'b', 'linewidth', 2);
scatter(nanmean(orientMatrix(:,strcmp(orientMatrixColIDs, 'TailX'))), nanmean(orientMatrix(:,strcmp(orientMatrixColIDs, 'TailY'))), 100, 'x', 'markeredgecolor', 'k', 'markerfacecolor', 'r', 'linewidth', 2);
title(sprintf('%s %s PreTrial = %i ms, PostTrial = %i ms', fileParts{1}, fileParts{2}, orientData.Params.PreTrialDuration*1000, orientData.Params.PostTrialDuration*1000));

trlStruct = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [0 0], 'PokeIn');
headSP = subplot(2,2,3);
tailSP = subplot(2,2,4);
linkaxes([scatAll, headSP, tailSP], 'xy');
for trl = 1:length(trlStruct)
    winStart = trlStruct(trl).PokeInIndex - (orientData.Params.PreTrialDuration*1000);
    winEnd = trlStruct(trl).PokeOutIndex + (orientData.Params.PostTrialDuration*1000);
    trlOrient = orientMatrix(winStart:winEnd,:);
    omx = trlOrient(sum(~isnan(trlOrient),2)==size(trlOrient,2),:);
    tempHead = plot(headSP, omx(:,strcmp(orientMatrixColIDs, 'HeadX')), omx(:,strcmp(orientMatrixColIDs, 'HeadY')), 'k');
    tempHead.Color(4) = 0.25;
    hold(headSP, 'on');
    tempTail = plot(tailSP, omx(:,strcmp(orientMatrixColIDs, 'TailX')), omx(:,strcmp(orientMatrixColIDs, 'TailY')), 'k');
    tempTail.Color(4) = 0.25;
    hold(tailSP, 'on');
    drawnow; 
    pause(1)
end

'