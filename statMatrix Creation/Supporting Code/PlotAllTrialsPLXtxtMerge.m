function PlotAllTrialsPLXtxtMerge(plxData)

trialDta = plxData.Raw;
txtDta = plxData.TextData;

figure
hold on;
corPerfScatX = [];
corPerfScatY = [];
inCorPerfScatX = [];
inCorPerfScatY = [];
plxHDscatX = [];
plxHDscatY = [];
for trl = 1:length(trialDta)
    curPokeMtrx = trialDta(trl).OdorPokesDurations;
%     curPokeMtrx = trialDta(trl).TrialPokesDurations;
    for p = 1:size(curPokeMtrx,1)
        if trialDta(trl).TranspositionDistance == 0
            line(curPokeMtrx(p,3:4), ones(1,2).*(length(trialDta)-trl+1), 'marker', '*', 'color', 'black');
        else
            line(curPokeMtrx(p,3:4), ones(1,2).*(length(trialDta)-trl+1), 'marker', '*', 'color', 'blue');
        end
    end
    if trialDta(trl).Performance == 1
        corPerfScatX = [corPerfScatX, txtDta(trl).PokeDur]; %#ok<*AGROW>
        corPerfScatY = [corPerfScatY, length(trialDta)-trl+1];
    else
        inCorPerfScatX = [inCorPerfScatX, txtDta(trl).PokeDur];
        inCorPerfScatY = [inCorPerfScatY, length(trialDta)-trl+1];
    end
    plxHDscatX = [plxHDscatX, trialDta(trl).PokeDuration];
    plxHDscatY = [plxHDscatY, length(trialDta)-trl+1];
end
corr = scatter(corPerfScatX, corPerfScatY,50,'r', 'linewidth', 2);
inCorr = scatter(inCorPerfScatX, inCorPerfScatY,100,'r', 'x', 'linewidth', 2);
scatter(plxHDscatX, plxHDscatY,20,'filled', 'g', 'o');
lgd = legend([corr inCorr], 'Correct', 'Incorrect');
title(lgd, 'Text File Performance');
line(repmat(plxData.Summary.TargetDuration, [1,2]), get(gca, 'ylim'));
set(gca, 'ytick', 1:length(trialDta), 'yticklabel', fliplr([trialDta.TrialNum]));
title({[plxData.Summary.PlxFile ' All trials']; 'Black lines = InSeq; Blue = OutSeq'; 'Red Circle = Text Timestamp, Green = Plexon'})
