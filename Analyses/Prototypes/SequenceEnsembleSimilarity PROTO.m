
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

load(fileNames{cellfun(@(a)~isempty(a), strfind(fileNames, 'BehaviorMatrix'))});
load(fileNames{cellfun(@(a)~isempty(a), strfind(fileNames, 'EnsembleMatrix'))});
%%
eventWindow = [-0.5 5];
pehBinSize = 0.125;

pokeInAlignedBehavMatrix = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, eventWindow, 'PokeIn');
pokeInEnsembles = ExtractTrialData_SM(pokeInAlignedBehavMatrix, ensembleMatrix(:,2:end));

pokeOutAlignedBehavMatrix = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, eventWindow, 'PokeOut');
pokeOutEnsembles = ExtractTrialData_SM(pokeOutAlignedBehavMatrix, ensembleMatrix(:,2:end));

binnedPokeIn = cell(size(pokeInEnsembles));
binnedPokeOut = cell(size(pokeOutEnsembles));
for trl = 1:length(pokeInEnsembles)
    [binnedPokeIn{trl}, pokeInBins] = RebinPEH2_SM(pokeInEnsembles{trl}, eventWindow, pehBinSize);
    [binnedPokeOut{trl}, pokeOutBins] = RebinPEH2_SM(pokeOutEnsembles{trl}, eventWindow, pehBinSize);
end

%%
trlSeqs = [pokeInAlignedBehavMatrix.SequenceNum];
seqs = unique(trlSeqs);
for seq = 1:length(seqs)
    curSeqTrls = pokeInAlignedBehavMatrix(trlSeqs==seq);
    if length(curSeqTrls) == 1
        continue
    else
        figure
        for trl1 = 1:length(curSeqTrls)
            trl1Num = curSeqTrls(trl1).TrialNum;
            trl1Nsmbl = binnedPokeIn{trl1Num};
            for trl2 = 1:length(curSeqTrls)
                trl2Num = curSeqTrls(trl2).TrialNum;
                trl2Nsmbl = binnedPokeIn{trl2Num};
                
                index = sub2ind([length(curSeqTrls), length(curSeqTrls)], trl1, trl2);
                subplot(length(curSeqTrls), length(curSeqTrls), index)
                imagesc(pokeInBins, pokeInBins, corr(trl1Nsmbl, trl2Nsmbl), [0 1]);
                hold on;
                line(get(gca, 'xlim'), [0 0], 'linewidth', 2, 'color', 'black');
                line(get(gca, 'xlim'), [curSeqTrls(trl1).PokeDuration, curSeqTrls(trl1).PokeDuration], 'linewidth', 2, 'color', 'black', 'linestyle', '--');
                line([0 0], get(gca, 'ylim'), 'linewidth', 2, 'color', 'black');
                line([curSeqTrls(trl2).PokeDuration, curSeqTrls(trl2).PokeDuration], get(gca, 'ylim'), 'linewidth', 2, 'color', 'black', 'linestyle', '--');
                ylabel(sprintf('%i - %i', curSeqTrls(trl1).Odor, curSeqTrls(trl1).Position));
                xlabel(sprintf('%i - %i', curSeqTrls(trl2).Odor, curSeqTrls(trl2).Position));
                if index==1
                    title(sprintf('Sequence %i', curSeqTrls(1).SequenceNum));
                end
            end
        end            
    end
    colormap jet
    drawnow
    pause(5)
end
