function PlotErrorTrialsPLXtxtMerge(plxData, txtData, plxSummary)

if nargin == 1 && isfield(plxData, 'Raw')
    trialDta = plxData.Raw;
    txtDta = plxData.TextData;
    targDur = plxData.Summary.TargetDuration;
    fileName = plxData.Summary.PlxFile;
elseif nargin == 3 && isfield(plxData, 'TrialNum')
    trialDta = plxData;
    txtDta = txtData;
    targDur = plxSummary.TargetDuration;
    fileName = plxSummary.PlxFile;
end

errTrlsLog = logical(cell2mat(arrayfun(@(a)~isempty(a.Errors), trialDta, 'uniformoutput', 0)));

errorTrls = trialDta(errTrlsLog);
errorTrlsTxt = txtDta(errTrlsLog);

figure; 
hold on;
corPerfScatX = [];
corPerfScatY = [];
inCorPerfScatX = [];
inCorPerfScatY = [];
plxHDscatX = [];
plxHDscatY = [];
for trl = 1:length(errorTrls)
    curPokeMtrx = errorTrls(trl).OdorPokesDurations;
    for p = 1:size(curPokeMtrx,1)
        if errorTrls(trl).TranspositionDistance == 0
            line(curPokeMtrx(p,3:4), ones(1,2).*(length(errorTrls)-trl+1), 'marker', '*', 'color', 'black');
        else
            line(curPokeMtrx(p,3:4), ones(1,2).*(length(errorTrls)-trl+1), 'marker', '*', 'color', 'blue');
        end
    end
    if errorTrls(trl).Performance == 1
        corPerfScatX = [corPerfScatX, errorTrlsTxt(trl).PokeDur]; %#ok<*AGROW>
        corPerfScatY = [corPerfScatY, length(errorTrls)-trl+1];
    else
        inCorPerfScatX = [inCorPerfScatX, errorTrlsTxt(trl).PokeDur];
        inCorPerfScatY = [inCorPerfScatY, length(errorTrls)-trl+1];
    end
    plxHDscatX = [plxHDscatX, errorTrls(trl).PokeDuration];
    plxHDscatY = [plxHDscatY, length(errorTrls)-trl+1];
end
corr = scatter(corPerfScatX, corPerfScatY,50,'r', 'linewidth', 2);
inCorr = scatter(inCorPerfScatX, inCorPerfScatY,100,'r', 'x', 'linewidth', 2);
scatter(plxHDscatX, plxHDscatY,20,'filled', 'g', 'o');
lgd = legend([corr inCorr], 'Correct', 'Incorrect');
title(lgd, 'Text File Performance');
line(repmat(targDur, [1,2]), get(gca, 'ylim'));
set(gca, 'ytick', 1:length(errorTrls), 'yticklabel', fliplr([errorTrls.TrialNum]));
title({[fileName ' "Error" trials']; 'Black lines = InSeq; Blue = OutSeq'; 'Red Circle = Text Timestamp, Green = Plexon'})
