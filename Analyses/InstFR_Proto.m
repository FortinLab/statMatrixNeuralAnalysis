
%%
origDir = cd;
[fileDir] = uigetdir(origDir);
if fileDir==0
    disp('Analysis Cancelled')
    return
else
    cd(fileDir)
end

files = dir(cd);
fileNames = {files.name};
matFiles = fileNames(cell2mat(cellfun(@(a)~isempty(a), strfind(fileNames, '.mat'), 'uniformoutput', 0)))';
statMatrixLog = false(size(matFiles));
for fl = 1:length(matFiles)
    variableInfo = who('-file', matFiles{fl});
    if sum(ismember(variableInfo, 'statMatrix'))==1
        statMatrixLog(fl) = true;
    end
end
smFiles = matFiles(statMatrixLog);

%% Load Relevant Data
load(fileNames{cellfun(@(a)~isempty(a), strfind(fileNames, 'BehaviorMatrix'))});
load(fileNames{cellfun(@(a)~isempty(a), strfind(fileNames, 'EnsembleMatrix'))});
numUni = length(ensembleUnitSummaries);

%% Define Variables
w = gausswin(125);
w = w/(length(w)*mode(diff(behavMatrix(:,1))));
%% Compute Instantaneous Firing Rate and Z-Normalize to Session
ensmblFR = nan(size(ensembleMatrix,1), size(ensembleMatrix,2)-1);
for uni = 2:size(ensembleMatrix,2)
    ensmblFR(:,uni-1) = zscore(conv(ensembleMatrix(:,uni), w, 'same'));
end

%% Extract and Analyze PokeIn Aligned Data
pokeInAlignedMatrix = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [-0.8, 0], 'PokeIn');
odorAlignedMatrix = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [0 0.5], 'Odor');
pokeOutAlignedMatrix = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [-0.5 0.2], 'PokeOut');
errorAlignedMatrix = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [-0.25 0.5], 'PokeOut');
rewardAlignedMatrix = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [-0.25 0.25], 'FrontReward');
isLog = [pokeInAlignedMatrix.TranspositionDistance]==0;
perfLog = [pokeInAlignedMatrix.Performance]==1;

piVect = -0.8:0.001:0;
odrVect = 0:0.001:0.5;
poVect = -0.5:0.001:0.2;
rwdVect = -0.25:0.001:0.25;
errVect = -0.25:0.001:0.5;

pokeInEnsmblFR = ExtractTrialData_SM(pokeInAlignedMatrix(isLog & perfLog), ensmblFR);
pokeInEnsmblRSTR = ExtractTrialData_SM(pokeInAlignedMatrix(isLog & perfLog), ensembleMatrix(:,2:end));
odorEnsmblFR = ExtractTrialData_SM(odorAlignedMatrix(isLog & perfLog), ensmblFR);
odorEnsmblRSTR = ExtractTrialData_SM(odorAlignedMatrix(isLog & perfLog), ensembleMatrix(:,2:end));
pokeOutEnsmblFR = ExtractTrialData_SM(pokeOutAlignedMatrix(isLog & perfLog), ensmblFR);
pokeOutEnsmblRSTR = ExtractTrialData_SM(pokeOutAlignedMatrix(isLog & perfLog), ensembleMatrix(:,2:end));
rewardEnsmblFR = ExtractTrialData_SM(rewardAlignedMatrix(isLog & perfLog), ensmblFR);
rewardEnsmblRSTR = ExtractTrialData_SM(rewardAlignedMatrix(isLog & perfLog), ensembleMatrix(:,2:end));
errEnsmblFR = ExtractTrialData_SM(errorAlignedMatrix(isLog & ~perfLog), ensmblFR);
errEnsmblRSTR = ExtractTrialData_SM(errorAlignedMatrix(isLog & ~perfLog), ensembleMatrix(:,2:end));
for uni = 1:numUni
    piTrlFR = cell2mat(cellfun(@(a){a(:,uni)}, pokeInEnsmblFR));
    piTrlRSTR = cell2mat(cellfun(@(a){a(:,uni)}, pokeInEnsmblRSTR));
    odrTrlFR = cell2mat(cellfun(@(a){a(:,uni)}, odorEnsmblFR));
    odrTrlRSTR = cell2mat(cellfun(@(a){a(:,uni)}, odorEnsmblRSTR));
    poTrlFR = cell2mat(cellfun(@(a){a(:,uni)}, pokeOutEnsmblFR));
    poTrlRSTR = cell2mat(cellfun(@(a){a(:,uni)}, pokeOutEnsmblRSTR));
    rwdTrlFR = cell2mat(cellfun(@(a){a(:,uni)}, rewardEnsmblFR));
    rwdTrlRSTR = cell2mat(cellfun(@(a){a(:,uni)}, rewardEnsmblRSTR));
    errTrlFR = cell2mat(cellfun(@(a){a(:,uni)}, errEnsmblFR));
    errTrlRSTR = cell2mat(cellfun(@(a){a(:,uni)}, errEnsmblRSTR));
    
    figure;
    pi = subplot(2,5,6);
    plot(piVect,piTrlFR);
    hold on;
    plot(piVect,mean(piTrlFR,2), 'linewidth', 2, 'color', 'black');
    title([ensembleMatrixColIDs{uni+1} ' Pre Trial']);
    
    piR = subplot(2,5,1);
    hold on;
    for trl = 1:size(piTrlRSTR,2)
        scatter(piVect(piTrlRSTR(:,trl)==1), ones(1,sum(piTrlRSTR(:,trl)==1)).*size(piTrlRSTR,2)+1-trl, '*k');
    end
    axis off;    
    linkaxes([pi,piR], 'x');
    axis(pi, 'tight');
        
    odr = subplot(2,5,7);
    plot(odrVect,odrTrlFR);
    hold on;
    plot(odrVect,mean(odrTrlFR,2), 'linewidth', 2,'color', 'black');
    title('Early Trial');
    
    odrR = subplot(2,5,2);
    hold on;
    for trl = 1:size(odrTrlRSTR,2)
        scatter(odrVect(odrTrlRSTR(:,trl)==1), ones(1,sum(odrTrlRSTR(:,trl)==1)).*size(odrTrlRSTR,2)+1-trl, '*k');
    end
    axis off;    
    linkaxes([odr,odrR], 'x');
    axis(odr, 'tight');
    
    po = subplot(2,5,8);
    plot(poVect, poTrlFR);
    hold on;
    plot(poVect,mean(poTrlFR,2), 'linewidth', 2, 'color', 'black');
    title('Late Trial');
    
    poR = subplot(2,5,3);
    hold on;
    for trl = 1:size(poTrlRSTR,2)
        scatter(poVect(poTrlRSTR(:,trl)==1), ones(1,sum(poTrlRSTR(:,trl)==1)).*size(poTrlRSTR,2)+1-trl, '*k');
    end
    axis off;    
    linkaxes([po,poR], 'x');
    axis(po, 'tight');
    
    rwd = subplot(2,5,9);
    plot(rwdVect,rwdTrlFR);
    hold on;
    plot(rwdVect,mean(rwdTrlFR,2), 'linewidth', 2, 'color', 'black');
    title('Reward');
    
    rwdR = subplot(2,5,4);
    hold on;
    for trl = 1:size(rwdTrlRSTR,2)
        scatter(rwdVect(rwdTrlRSTR(:,trl)==1), ones(1,sum(rwdTrlRSTR(:,trl)==1)).*size(rwdTrlRSTR,2)+1-trl, '*k');
    end
    axis off;    
    linkaxes([rwd,rwdR], 'x');
    axis(rwd, 'tight');
    
    err = subplot(2,5,10);
    plot(errVect,errTrlFR);
    hold on;
    plot(errVect,mean(errTrlFR,2), 'linewidth', 2, 'color', 'black');
    title('Error');
    
    errR = subplot(2,5,5);
    hold on;
    for trl = 1:size(errTrlRSTR,2)
        scatter(errVect(errTrlRSTR(:,trl)==1), ones(1,sum(errTrlRSTR(:,trl)==1)).*size(errTrlRSTR,2)+1-trl, '*k');
    end
    axis off;    
    linkaxes([err,errR], 'x');
    axis(err, 'tight');
    
%     set([pi,odr,po,rwd,err], 'yscale', 'log');
    linkaxes([pi,odr,po,rwd,err], 'y');
    drawnow
end

%%


pokeOutAlignedMatrix = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [-0.5, 0.5], 'PokeOut');
rewardAlignedMatrix = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [-0.5 0.5], 'FrontReward');





% Extract time bins
eventTimes = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [0 0], 'PokeIn');
trialTimeBins = ExtractTrialData_SM(pokeInAlignedMatrix, behavMatrix(:,1));
eventTimes = ExtractTrialData_SM(eventTimes, behavMatrix(:,1));
eventTimeBins = cellfun(@(a,b) a-b, trialTimeBins, eventTimes, 'uniformoutput',0);
frstNonMptTrl = find(cellfun(@(a)~isempty(a), eventTimeBins),1, 'first');
eventTimeBins = eventTimeBins{frstNonMptTrl};
% Extract Data
pokeOutEnsmblFR = ExtractTrialData_SM(pokeOutAlignedMatrix, ensmblFR);

figure
pi = subplot(2,1,1);
plot(eventTimeBins,mean(cell2mat(reshape(pokeInEnsmblFR, [1,1,length(pokeInEnsmblFR)])),3));
po = subplot(2,1,2);
plot(eventTimeBins, mean(cell2mat(reshape(pokeOutEnsmblFR, [1,1,length(pokeOutEnsmblFR)])),3));
linkaxes([pi, po], 'y');

%% Identify and load LFP data
curFile = fileNames{cellfun(@(a)~isempty(a),regexp(fileNames, sprintf('([A-Z]*[0-9]*)_([A-Z | a-z]*)([0-9]*)_([A-Z | a-z]*)([0-9]*)_%s_SM.mat', curTet)))};
load(curFile);
hilbColLog = cellfun(@(a)~isempty(a),strfind(statMatrixColIDs, '_HilbVals'));
hilbCols = find(hilbColLog);
hilbCols = hilbCols(2:end);
bandNames = cellfun(@(b)b{3}, cellfun(@(a)strsplit(a, '_'), statMatrixColIDs(hilbColLog), 'uniformoutput',0), 'uniformoutput',0);
bands = bandNames(2:end);
