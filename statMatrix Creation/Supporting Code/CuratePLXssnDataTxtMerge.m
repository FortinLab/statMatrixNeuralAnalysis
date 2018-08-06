function [plxData] = CuratePLXssnDataTxtMerge(plxData)
%%
%% Check inputs.... tbd
if nargin==0
    uigetfile(blah)
end
%% Extract data
plxSession = plxData.Raw;
plxText = plxData.TextData;
identifier = plxData.Summary.Identifier;

errTrials = find(logical(cell2mat(arrayfun(@(a)~isempty(a.Errors), plxSession, 'uniformoutput',0))));

outfile = fopen(sprintf('%s_CuratePLX_%i.txt', plxData.Summary.PlxFile(1:end-4), identifier), 'wt');
fprintf(outfile, '%s txtDta file merger PLX timestamp curation summary\n', plxData.Summary.PlxFile);


%% Run through the error trials
flRmvlLog = false(size(plxSession));
for trl = 1:length(errTrials)
    h = figure('position', [680 750 560 210]);
    hold on
    curPokeMtrx = plxSession(errTrials(trl)).OdorPokesDurations;
%     curPokeMtrx = plxSession(errTrials(trl)).TrialPokesDurations;
    for p = 1:size(curPokeMtrx,1)
        if plxSession(errTrials(trl)).TranspositionDistance == 0
            line(curPokeMtrx(p,3:4), ones(1,2), 'marker', '*', 'color', 'black');
        else
            line(curPokeMtrx(p,3:4), ones(1,2), 'marker', '*', 'color', 'blue');
        end
    end
    if plxSession(errTrials(trl)).Performance == 1
        scatter(plxText(errTrials(trl)).PokeDur, 1,50,'r', 'linewidth', 2);
    else
        scatter(plxText(errTrials(trl)).PokeDur, 1,100,'r', 'x', 'linewidth', 2);
    end
    scatter(plxSession(errTrials(trl)).PokeDuration, 1, 20,'filled', 'g', 'o');
    set(gca, 'ylim', [0 2], 'xlim', [0 3]);
    line(repmat(plxData.Summary.TargetDuration, [1,2]), get(gca, 'ylim'));
    
    quest = questdlg(sprintf('How to handle trial #%i?', plxSession(errTrials(trl)).TrialNum),...
        sprintf('Review Trial #%i', plxSession(errTrials(trl)).TrialNum),...
        'Use .PLX timestamp',...
        'Use .TXT timestamp',...
        'Remove Trial',...
        'Use .PLX timestamp');
    close(h);
    if isempty(quest)
        disp('Curation aborted');
        fprintf(outfile, 'Curation aborted examining trial #%i (%i/%i)\n', plxSession(errTrials(trl)).TrialNum,...
            trl, length(errTrials));
        return
    end
    switch quest
        case 'Use .PLX timestamp'
            fprintf(outfile, 'PLX timestamp used for trial #%i\n', plxSession(errTrials(trl)).TrialNum);
                        
        case 'Use .TXT timestamp'
            fprintf(outfile, 'TXT timestamp used for trial #%i\n', plxSession(errTrials(trl)).TrialNum);
            pokeDiffMtx = curPokeMtrx(:,2) - curPokeMtrx(:,1)';
            if sum(sum(round(plxText(errTrials(trl)).PokeDur,4) == round(pokeDiffMtx,4)))==1
                [withdrawNum,~] = find(round(plxText(errTrials(trl)).PokeDur,4) == round(pokeDiffMtx,4));
            elseif sum(sum(round(plxText(errTrials(trl)).PokeDur,5) == round(pokeDiffMtx,5)))==1
                [withdrawNum,~] = find(round(plxText(errTrials(trl)).PokeDur,5) == round(pokeDiffMtx,5));
            else
                error('Cannot find .txt poke duration... something is wrong');
            end
            plxSession(errTrials(trl)).OdorPokeWithdrawTime = curPokeMtrx(withdrawNum,2);
            plxSession(errTrials(trl)).PokeDuration = plxSession(errTrials(trl)).OdorPokeWithdrawTime - plxSession(errTrials(trl)).OdorTrigPokeTime;
        case 'Remove Trial'
            fprintf(outfile, 'REMOVED trial #%i\n', plxSession(errTrials(trl)).TrialNum);
            flRmvlLog(errTrials(trl)) = true;
    end    
end
plxSession(flRmvlLog) = [];
plxText(flRmvlLog) = [];
PlotErrorTrialsPLXtxtMerge(plxSession, plxText, plxData.Summary)

%%
plxData.CuratedPLX = plxSession;
plxData.CuratedTXT = plxText;
fclose(outfile);
uisave('plxData', [plxData.Summary.PlxFile(1:end-4) '_plxData.mat']);