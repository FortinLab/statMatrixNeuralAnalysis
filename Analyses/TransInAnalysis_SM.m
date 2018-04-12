function TransInAnalysis_SM
%% IN DEVELOPMENT %%


%%
close all
clear all
%%
pehBinSize = 0.25;
eventWindow = [-1 1];
%%
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

%% Load Relevant Data
load(fileNames{cellfun(@(a)~isempty(a), strfind(fileNames, 'BehaviorMatrix'))});
load(fileNames{cellfun(@(a)~isempty(a), strfind(fileNames, 'EnsembleMatrix'))});

%%
pokeInAlignedBehavMatrix = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, eventWindow, 'PokeIn');

outSeqLog = ~([pokeInAlignedBehavMatrix.TranspositionDistance]==0);
trialPerfLog = [pokeInAlignedBehavMatrix.Performance]==1;
fourthPosLog = [pokeInAlignedBehavMatrix.Position]==pokeInAlignedBehavMatrix(1).SeqLength;

outSeqCorrNonFourthPosTrls = find(outSeqLog&trialPerfLog&~(fourthPosLog));
outSeqCorrOdor = [pokeInAlignedBehavMatrix(outSeqCorrNonFourthPosTrls).Odor];
transInTrls = outSeqCorrNonFourthPosTrls+1;
% Remove those trials they got wrong
outSeqCorrOdor([pokeInAlignedBehavMatrix(transInTrls).Performance]==0) = [];
transInTrls([pokeInAlignedBehavMatrix(transInTrls).Performance]==0) = [];
transInTrlsPos = [pokeInAlignedBehavMatrix(transInTrls).Position];

goodTrials = transInTrls;
odorIDs = outSeqCorrOdor;
posIDs = transInTrlsPos;

%%
pokeInAlignedBehavMatrix = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, eventWindow, 'PokeIn');

trialPerfLog = [pokeInAlignedBehavMatrix.Performance]==0;
firstPosLog = [pokeInAlignedBehavMatrix.Position]==1;


goodTrials = find(trialPerfLog&~firstPosLog);
odorIDs = [pokeInAlignedBehavMatrix(goodTrials).Odor];
posIDs = [pokeInAlignedBehavMatrix(goodTrials).Position];

%%
fRatDiff = cell(length(ensembleMatrixColIDs)-1,1); %#ok<USENS>
for u = 2:length(ensembleMatrixColIDs) 
    curUni = ensembleMatrixColIDs{u};
    curUniSpkData = ensembleMatrix(:,u); %#ok<NODEF>
    curUniData = ExtractTrialData_SM(pokeInAlignedBehavMatrix(goodTrials), curUniSpkData);
    noSpkLog = cellfun(@(a)isempty(a), curUniData);
    curUniData(noSpkLog) = [];
    curUniPEH = cell(length(curUniData),1);
    if ~isempty(curUniData)
        for trl = 1:length(curUniPEH)
            [curUniPEH{trl,1}, ~] = RebinPEH_SM(curUniData{trl}, eventWindow, pehBinSize);
        end
        curUniUnpacked = cell2mat(curUniPEH);
        odorF = nan(1,size(curUniUnpacked,2));
        positionF = nan(1,size(curUniUnpacked,2));
        for ndx = 1:size(curUniUnpacked,2)
            [~,odrTable,~] = anova1(curUniUnpacked(:,ndx), odorIDs', 'off');
            curOdorF = odrTable{2,5}; 
            if ~isempty(curOdorF) && curOdorF>0 && curOdorF<5
                odorF(ndx) = curOdorF;
            elseif curOdorF>5
                odorF(ndx) = 5;
            else
                odorF(ndx) = 0;
            end
            [~,posTable,~] = anova1(curUniUnpacked(:,ndx), posIDs', 'off');
            curPosF = posTable{2,5};
            if ~isempty(curPosF) && curPosF>0 && curPosF<5
                positionF(ndx) = curPosF;
            elseif curPosF>5
                positionF(ndx) = 5;
            else
                positionF(ndx) = 0;
            end
        end
    end
    fRatDiff{u-1} = positionF-odorF;
end

fRatDiffMtx = cell2mat(fRatDiff);
transInMean = nanmean(fRatDiffMtx);
transInSEM = nanstd(fRatDiffMtx,0,1)./(sum(~isnan(fRatDiffMtx))-1);
xTicks = (eventWindow(1):pehBinSize:(eventWindow(2)-pehBinSize))+(pehBinSize/2);

figure
PlotLineAndFilledError(xTicks, transInMean, transInSEM, 'blue');

