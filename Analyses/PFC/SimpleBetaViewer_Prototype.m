clc
clear all

%%
piWindow = [-0.5 1.5];
poWindow = [-1.5 1.5];
%%
smPath = uigetdir;
cd(smPath);
files = dir(smPath);
fileNames = {files.name};
% Load the behavior matrix file for poke events plots
behMatFile = fileNames{cellfun(@(a)~isempty(a), strfind(fileNames, 'BehaviorMatrix'))};
% nsmblMatFile = fileNames{cellfun(@(a)~isempty(a), strfind(fileNames, 'EnsembleMatrix'))};
load([smPath '\' behMatFile]);
% load([smPath '\' nsmblMatFile]);
% Identify list of all statMatrix files
smFileList = fileNames(cellfun(@(a)~isempty(a), regexp(fileNames, '_SM\>')))';
% Organize Trial Data and Define Trial Logicals
tdPI = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, piWindow, 'PokeIn');
tdPO = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, poWindow, 'PokeOut');
tdRS = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, poWindow, 'RewardSignal');
isLog = [tdPI.TranspositionDistance]==0;
perfLog = [tdPI.Performance]==1;
p1Log = [tdPI.Position]==1;
p2Log = [tdPI.Position]==2;
p3Log = [tdPI.Position]==3;
p4Log = [tdPI.Position]==4;
pokeDur = [tdPI.PokeDuration];
wdLat = [tdPI.WithdrawLatency];

targDur = nanmean([tdPI.RewardSignalIndex]-[tdPI.PokeInIndex])/1000;
targVar = nanstd([tdPI.RewardSignalIndex]-[tdPI.PokeInIndex])/1000;

% Define time period of the session based on first and last trial to remove
% sleep periods from the z-score calculation
ssnLog = false(size(behavMatrix,1),1);
ssnLog(tdPI(1).PokeInIndex-500:max([tdPI.PokeOutIndex])+500) = true;

%%
wlCorrVals = cell(size(smFileList));
tetIDs = cell(size(smFileList));
for f = 1:length(smFileList)
    % Load File
    load([smPath '\' smFileList{f}]);
    % Identify Beta, get envelope and z-score
    betaCol = cellfun(@(a)~isempty(a), regexp(statMatrixColIDs, '_Beta\>'));
    tet = strsplit(statMatrixColIDs{betaCol}, '_');
    tetIDs(f) = tet(1);
%     env = abs(hilbert(statMatrix(:,betaCol)));
    env = abs(hilbert(statMatrix(:,betaCol))).^2;
    zEnv = zeros(size(statMatrix,1),1);
    zEnv(ssnLog) = zscore(env(ssnLog));
    
    % InSeq v OutSeq Beta Power Overall
    iscBetaPI = cell2mat(ExtractTrialData_SM(tdPI(isLog & perfLog), zEnv));
    oscBetaPI = cell2mat(ExtractTrialData_SM(tdPI(~isLog & perfLog), zEnv));
    iscBetaPO = cell2mat(ExtractTrialData_SM(tdPO(isLog & perfLog), zEnv));
    oscBetaPO = cell2mat(ExtractTrialData_SM(tdPO(~isLog & perfLog), zEnv));
%     figure;
%     piSP = subplot(6,1,1:2);
%     line(piWindow, [0 0], 'linestyle', ':', 'color', 'k');
%     hold on;
%     plot(piWindow(1):1/1000:piWindow(2), mean(iscBetaPI,2), 'color', 'blue');
%     patch('YData', [mean(iscBetaPI,2)+std(iscBetaPI,1,2)/sqrt(size(iscBetaPI,2)-1); flipud(mean(iscBetaPI,2)-std(iscBetaPI,1,2)/sqrt(size(iscBetaPI,2)-1))],...
%         'XData', [piWindow(1):1/1000:piWindow(2), piWindow(2):-1/1000:piWindow(1)],...
%         'FaceColor', 'blue', 'FaceAlpha', .3, 'EdgeColor', 'blue', 'EdgeAlpha', 0.5);
%     plot(piWindow(1):1/1000:piWindow(2), mean(oscBetaPI,2), 'color', 'red');
%     patch('YData', [mean(oscBetaPI,2)+std(oscBetaPI,1,2)/sqrt(size(oscBetaPI,2)-1); flipud(mean(oscBetaPI,2)-std(oscBetaPI,1,2)/sqrt(size(oscBetaPI,2)-1))],...
%         'XData', [piWindow(1):1/1000:piWindow(2), piWindow(2):-1/1000:piWindow(1)],...
%         'FaceColor', 'red', 'FaceAlpha', .3, 'EdgeColor', 'red', 'EdgeAlpha', 0.5);
%     title('Poke In Aligned');
%     subplot(6,1,3)
%     histogram(pokeDur(isLog & perfLog),piWindow(1):0.01:piWindow(2), 'normalization', 'count');
%     hold on;
%     histogram(pokeDur(~isLog & perfLog),piWindow(1):0.01:piWindow(2), 'normalization', 'count');
%     axis tight
%     set(gca, 'xlim', get(piSP, 'xlim'));
%     poSP = subplot(6,1,4:5);
%     line(poWindow, [0 0], 'linestyle', ':', 'color', 'k');
%     hold on;
%     plot(poWindow(1):1/1000:poWindow(2), mean(iscBetaPO,2), 'color', 'blue');
%     patch('YData', [mean(iscBetaPO,2)+std(iscBetaPO,1,2)/sqrt(size(iscBetaPO,2)-1); flipud(mean(iscBetaPO,2)-std(iscBetaPO,1,2)/sqrt(size(iscBetaPO,2)-1))],...
%         'XData', [poWindow(1):1/1000:poWindow(2), poWindow(2):-1/1000:poWindow(1)],...
%         'FaceColor', 'blue', 'FaceAlpha', .3, 'EdgeColor', 'blue', 'EdgeAlpha', 0.5);
%     plot(poWindow(1):1/1000:poWindow(2), mean(oscBetaPO,2), 'color', 'red');
%     patch('YData', [mean(oscBetaPO,2)+std(oscBetaPO,1,2)/sqrt(size(oscBetaPO,2)-1); flipud(mean(oscBetaPO,2)-std(oscBetaPO,1,2)/sqrt(size(oscBetaPO,2)-1))],...
%         'XData', [poWindow(1):1/1000:poWindow(2), poWindow(2):-1/1000:poWindow(1)],...
%         'FaceColor', 'red', 'FaceAlpha', .3, 'EdgeColor', 'red', 'EdgeAlpha', 0.5);
%     title('Poke Out Aligned');
%     subplot(6,1,6)
%     histogram(wdLat(isLog & perfLog)*-1, poWindow(1):0.01:poWindow(2));
%     hold on;
%     histogram(pokeDur(~isLog & perfLog)*-1, poWindow(1):0.01:poWindow(2));
%     axis tight
%     set(gca, 'xlim', get(poSP, 'xlim'));
%     annotation('textbox', 'position', [0.025 0.025 0.7 0.05], 'String', smFileList{f},...
%         'linestyle', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
%     drawnow;
%     close gcf
%     
%     % InSeq & OutSeq Beta Broken down by Position
%     figure;
%     sp(1) = subplot(2,4,1);
%     iscBetaPI1 = cell2mat(ExtractTrialData_SM(tdPI(isLog & perfLog & p1Log), zEnv));
%     line(piWindow, [0 0], 'linestyle', ':', 'color', 'k');
%     hold on;
%     plot(piWindow(1):1/1000:piWindow(2), mean(iscBetaPI1,2), 'color', 'blue');
%     patch('YData', [mean(iscBetaPI1,2)+std(iscBetaPI1,1,2)/sqrt(size(iscBetaPI1,2)-1); flipud(mean(iscBetaPI1,2)-std(iscBetaPI1,1,2)/sqrt(size(iscBetaPI1,2)-1))],...
%         'XData', [piWindow(1):1/1000:piWindow(2), piWindow(2):-1/1000:piWindow(1)],...
%         'FaceColor', 'blue', 'FaceAlpha', .3, 'EdgeColor', 'blue', 'EdgeAlpha', 0.5);
%     ylabel('Poke In Aligned');
%     title('Position 1');
%     sp(2) = subplot(2,4,2);
%     iscBetaPI2 = cell2mat(ExtractTrialData_SM(tdPI(isLog & perfLog & p2Log), zEnv));
%     oscBetaPI2 = cell2mat(ExtractTrialData_SM(tdPI(~isLog & perfLog & p2Log), zEnv));
%     line(piWindow, [0 0], 'linestyle', ':', 'color', 'k');
%     hold on;
%     plot(piWindow(1):1/1000:piWindow(2), mean(iscBetaPI2,2), 'color', 'blue');
%     patch('YData', [mean(iscBetaPI2,2)+std(iscBetaPI2,1,2)/sqrt(size(iscBetaPI2,2)-1); flipud(mean(iscBetaPI2,2)-std(iscBetaPI2,1,2)/sqrt(size(iscBetaPI2,2)-1))],...
%         'XData', [piWindow(1):1/1000:piWindow(2), piWindow(2):-1/1000:piWindow(1)],...
%         'FaceColor', 'blue', 'FaceAlpha', .3, 'EdgeColor', 'blue', 'EdgeAlpha', 0.5);
%     plot(piWindow(1):1/1000:piWindow(2), mean(oscBetaPI2,2), 'color', 'red');
%     patch('YData', [mean(oscBetaPI2,2)+std(oscBetaPI2,1,2)/sqrt(size(oscBetaPI2,2)-1); flipud(mean(oscBetaPI2,2)-std(oscBetaPI2,1,2)/sqrt(size(oscBetaPI2,2)-1))],...
%         'XData', [piWindow(1):1/1000:piWindow(2), piWindow(2):-1/1000:piWindow(1)],...
%         'FaceColor', 'red', 'FaceAlpha', .3, 'EdgeColor', 'red', 'EdgeAlpha', 0.5);
%     title('Position 2');
%     sp(3) = subplot(2,4,3);
%     iscBetaPI3 = cell2mat(ExtractTrialData_SM(tdPI(isLog & perfLog & p3Log), zEnv));
%     oscBetaPI3 = cell2mat(ExtractTrialData_SM(tdPI(~isLog & perfLog & p3Log), zEnv));
%     line(piWindow, [0 0], 'linestyle', ':', 'color', 'k');
%     hold on;
%     plot(piWindow(1):1/1000:piWindow(2), mean(iscBetaPI3,2), 'color', 'blue');
%     patch('YData', [mean(iscBetaPI3,2)+std(iscBetaPI3,1,2)/sqrt(size(iscBetaPI3,2)-1); flipud(mean(iscBetaPI3,2)-std(iscBetaPI3,1,2)/sqrt(size(iscBetaPI3,2)-1))],...
%         'XData', [piWindow(1):1/1000:piWindow(2), piWindow(2):-1/1000:piWindow(1)],...
%         'FaceColor', 'blue', 'FaceAlpha', .3, 'EdgeColor', 'blue', 'EdgeAlpha', 0.5);
%     plot(piWindow(1):1/1000:piWindow(2), mean(oscBetaPI3,2), 'color', 'red');
%     patch('YData', [mean(oscBetaPI3,2)+std(oscBetaPI3,1,2)/sqrt(size(oscBetaPI3,2)-1); flipud(mean(oscBetaPI3,2)-std(oscBetaPI3,1,2)/sqrt(size(oscBetaPI3,2)-1))],...
%         'XData', [piWindow(1):1/1000:piWindow(2), piWindow(2):-1/1000:piWindow(1)],...
%         'FaceColor', 'red', 'FaceAlpha', .3, 'EdgeColor', 'red', 'EdgeAlpha', 0.5);
%     title('Position 3');
%     sp(4) = subplot(2,4,4);
%     iscBetaPI4 = cell2mat(ExtractTrialData_SM(tdPI(isLog & perfLog & p4Log), zEnv));
%     oscBetaPI4 = cell2mat(ExtractTrialData_SM(tdPI(~isLog & perfLog & p4Log), zEnv));
%     line(piWindow, [0 0], 'linestyle', ':', 'color', 'k');
%     hold on;
%     plot(piWindow(1):1/1000:piWindow(2), mean(iscBetaPI4,2), 'color', 'blue');
%     patch('YData', [mean(iscBetaPI4,2)+std(iscBetaPI4,1,2)/sqrt(size(iscBetaPI4,2)-1); flipud(mean(iscBetaPI4,2)-std(iscBetaPI4,1,2)/sqrt(size(iscBetaPI4,2)-1))],...
%         'XData', [piWindow(1):1/1000:piWindow(2), piWindow(2):-1/1000:piWindow(1)],...
%         'FaceColor', 'blue', 'FaceAlpha', .3, 'EdgeColor', 'blue', 'EdgeAlpha', 0.5);
%     plot(piWindow(1):1/1000:piWindow(2), mean(oscBetaPI4,2), 'color', 'red');
%     patch('YData', [mean(oscBetaPI4,2)+std(oscBetaPI4,1,2)/sqrt(size(oscBetaPI4,2)-1); flipud(mean(oscBetaPI4,2)-std(oscBetaPI4,1,2)/sqrt(size(oscBetaPI4,2)-1))],...
%         'XData', [piWindow(1):1/1000:piWindow(2), piWindow(2):-1/1000:piWindow(1)],...
%         'FaceColor', 'red', 'FaceAlpha', .3, 'EdgeColor', 'red', 'EdgeAlpha', 0.5);
%     title('Position 4');
%     sp(5) = subplot(2,4,5);
%     iscBetaPO1 = cell2mat(ExtractTrialData_SM(tdPO(isLog & perfLog & p1Log), zEnv));
%     line(poWindow, [0 0], 'linestyle', ':', 'color', 'k');
%     hold on;
%     plot(poWindow(1):1/1000:poWindow(2), mean(iscBetaPO1,2), 'color', 'blue');
%     patch('YData', [mean(iscBetaPO1,2)+std(iscBetaPO1,1,2)/sqrt(size(iscBetaPO1,2)-1); flipud(mean(iscBetaPO1,2)-std(iscBetaPO1,1,2)/sqrt(size(iscBetaPO1,2)-1))],...
%         'XData', [poWindow(1):1/1000:poWindow(2), poWindow(2):-1/1000:poWindow(1)],...
%         'FaceColor', 'blue', 'FaceAlpha', .3, 'EdgeColor', 'blue', 'EdgeAlpha', 0.5);
%     ylabel('Poke Out Aligned');
%     sp(6) = subplot(2,4,6);
%     iscBetaPO2 = cell2mat(ExtractTrialData_SM(tdPO(isLog & perfLog & p2Log), zEnv));
%     oscBetaPO2 = cell2mat(ExtractTrialData_SM(tdPO(~isLog & perfLog & p2Log), zEnv));
%     line(poWindow, [0 0], 'linestyle', ':', 'color', 'k');
%     hold on;
%     plot(poWindow(1):1/1000:poWindow(2), mean(iscBetaPO2,2), 'color', 'blue');
%     patch('YData', [mean(iscBetaPO2,2)+std(iscBetaPO2,1,2)/sqrt(size(iscBetaPO2,2)-1); flipud(mean(iscBetaPO2,2)-std(iscBetaPO2,1,2)/sqrt(size(iscBetaPO2,2)-1))],...
%         'XData', [poWindow(1):1/1000:poWindow(2), poWindow(2):-1/1000:poWindow(1)],...
%         'FaceColor', 'blue', 'FaceAlpha', .3, 'EdgeColor', 'blue', 'EdgeAlpha', 0.5);
%     plot(poWindow(1):1/1000:poWindow(2), mean(oscBetaPO2,2), 'color', 'red');
%     patch('YData', [mean(oscBetaPO2,2)+std(oscBetaPO2,1,2)/sqrt(size(oscBetaPO2,2)-1); flipud(mean(oscBetaPO2,2)-std(oscBetaPO2,1,2)/sqrt(size(oscBetaPO2,2)-1))],...
%         'XData', [poWindow(1):1/1000:poWindow(2), poWindow(2):-1/1000:poWindow(1)],...
%         'FaceColor', 'red', 'FaceAlpha', .3, 'EdgeColor', 'red', 'EdgeAlpha', 0.5);
%     sp(7) = subplot(2,4,7);
%     iscBetaPO3 = cell2mat(ExtractTrialData_SM(tdPO(isLog & perfLog & p3Log), zEnv));
%     oscBetaPO3 = cell2mat(ExtractTrialData_SM(tdPO(~isLog & perfLog & p3Log), zEnv));
%     line(poWindow, [0 0], 'linestyle', ':', 'color', 'k');
%     hold on;
%     plot(poWindow(1):1/1000:poWindow(2), mean(iscBetaPO3,2), 'color', 'blue');
%     patch('YData', [mean(iscBetaPO3,2)+std(iscBetaPO3,1,2)/sqrt(size(iscBetaPO3,2)-1); flipud(mean(iscBetaPO3,2)-std(iscBetaPO3,1,2)/sqrt(size(iscBetaPO3,2)-1))],...
%         'XData', [poWindow(1):1/1000:poWindow(2), poWindow(2):-1/1000:poWindow(1)],...
%         'FaceColor', 'blue', 'FaceAlpha', .3, 'EdgeColor', 'blue', 'EdgeAlpha', 0.5);
%     plot(poWindow(1):1/1000:poWindow(2), mean(oscBetaPO3,2), 'color', 'red');
%     patch('YData', [mean(oscBetaPO3,2)+std(oscBetaPO3,1,2)/sqrt(size(oscBetaPO3,2)-1); flipud(mean(oscBetaPO3,2)-std(oscBetaPO3,1,2)/sqrt(size(oscBetaPO3,2)-1))],...
%         'XData', [poWindow(1):1/1000:poWindow(2), poWindow(2):-1/1000:poWindow(1)],...
%         'FaceColor', 'red', 'FaceAlpha', .3, 'EdgeColor', 'red', 'EdgeAlpha', 0.5);
%     sp(8) = subplot(2,4,8);
%     iscBetaPO4 = cell2mat(ExtractTrialData_SM(tdPO(isLog & perfLog & p4Log), zEnv));
%     oscBetaPO4 = cell2mat(ExtractTrialData_SM(tdPO(~isLog & perfLog & p4Log), zEnv));
%     line(poWindow, [0 0], 'linestyle', ':', 'color', 'k');
%     hold on;
%     plot(poWindow(1):1/1000:poWindow(2), mean(iscBetaPO4,2), 'color', 'blue');
%     patch('YData', [mean(iscBetaPO4,2)+std(iscBetaPO4,1,2)/sqrt(size(iscBetaPO4,2)-1); flipud(mean(iscBetaPO4,2)-std(iscBetaPO4,1,2)/sqrt(size(iscBetaPO4,2)-1))],...
%         'XData', [poWindow(1):1/1000:poWindow(2), poWindow(2):-1/1000:poWindow(1)],...
%         'FaceColor', 'blue', 'FaceAlpha', .3, 'EdgeColor', 'blue', 'EdgeAlpha', 0.5);
%     plot(poWindow(1):1/1000:poWindow(2), mean(oscBetaPO4,2), 'color', 'red');
%     patch('YData', [mean(oscBetaPO4,2)+std(oscBetaPO4,1,2)/sqrt(size(oscBetaPO4,2)-1); flipud(mean(oscBetaPO4,2)-std(oscBetaPO4,1,2)/sqrt(size(oscBetaPO4,2)-1))],...
%         'XData', [poWindow(1):1/1000:poWindow(2), poWindow(2):-1/1000:poWindow(1)],...
%         'FaceColor', 'red', 'FaceAlpha', .3, 'EdgeColor', 'red', 'EdgeAlpha', 0.5);
%     linkaxes(sp(1:4), 'x');
%     linkaxes(sp(5:8), 'x');
%     set(sp(1), 'xlim', piWindow);
%     set(sp(8), 'xlim', poWindow);
%     linkaxes(sp, 'y');
%     annotation('textbox', 'position', [0.025 0.025 0.7 0.05], 'String', smFileList{f},...
%         'linestyle', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
%     drawnow;

    % Examine Relationship between Withdrawal Latency and Beta Power
    iscBetaRS = cell2mat(ExtractTrialData_SM(tdRS(isLog & perfLog & ~p1Log), zEnv));
%     iscBetaRS = cell2mat(ExtractTrialData_SM(tdPO(isLog & perfLog), zEnv));
    pdTrls = pokeDur(isLog & perfLog & ~p1Log);
    wlTrls = wdLat(isLog & perfLog & ~p1Log);
    figure;
    
    % Plot Hold Durations
    subplot(2,4,[1,5]);
    scatter(ones(1,length(pdTrls))+normrnd(0,0.1, [1,length(pdTrls)]), pdTrls, 'marker', '*', 'markeredgecolor', 'k');
    hold on;
    line(get(gca, 'xlim'), [targDur targDur], 'color', 'k', 'linewidth', 2);
    line(get(gca, 'xlim'), [targDur+targVar targDur+targVar], 'color', 'k', 'linestyle', ':');
    line(get(gca, 'xlim'), [targDur-targVar targDur-targVar], 'color', 'k', 'linestyle', ':');
    set(gca, 'ylim', [0 max(pdTrls)+0.2]);
    title ('Poke Dur');
    
    % Plot Reaction Times
    subplot(2,4,[2,6]);
    scatter(ones(1,length(wlTrls))+normrnd(0,0.1, [1,length(wlTrls)]), wlTrls, 'marker', '*', 'markeredgecolor', 'k');
    hold on;
    line(get(gca, 'xlim'), [0 0], 'color', 'k', 'linewidth', 2);
    title('WD Lat');
    
    % Plot Beta Power relative to Reward Signal
    subplot(2,4,3:4)
    line(poWindow, [0 0], 'linestyle', ':', 'color', 'k');
    hold on;
    plot(poWindow(1):1/1000:poWindow(2), mean(iscBetaRS,2), 'color', 'blue');
    patch('YData', [mean(iscBetaRS,2)+std(iscBetaRS,1,2)/sqrt(size(iscBetaRS,2)-1); flipud(mean(iscBetaRS,2)-std(iscBetaRS,1,2)/sqrt(size(iscBetaRS,2)-1))],...
        'XData', [poWindow(1):1/1000:poWindow(2), poWindow(2):-1/1000:poWindow(1)],...
        'FaceColor', 'blue', 'FaceAlpha', .3, 'EdgeColor', 'blue', 'EdgeAlpha', 0.5);
   
    % Plot Correlation between Beta Power and behavior (overall Poke Duration and Reward Latency) over time
    pdCorr = nan(size(iscBetaRS,1),1);
    wlCorr = nan(size(iscBetaRS,1),1);
    for t = 1:size(iscBetaRS,1)
%         pdCorr(t) = corr(iscBetaRS(t,:)', pdTrls');
        pdCorr(t) = corr(iscBetaRS(t,:)', pdTrls', 'Type', 'Spearman');
%         wlCorr(t) = corr(iscBetaRS(t,:)', wlTrls');
        wlCorr(t) = corr(iscBetaRS(t,:)', wlTrls', 'Type', 'Spearman');
    end
    wlCorrVals{f} = wlCorr;
    subplot(2,4,7:8);
    plot(poWindow(1):1/1000:poWindow(2), pdCorr, 'color', 'k');
    hold on; 
    plot(poWindow(1):1/1000:poWindow(2), wlCorr, 'color', 'k', 'linestyle', '-.');
    legend('Poke Duration', 'Withdrawal Latency', 'location', 'southoutside');
    line(poWindow, [0 0], 'color', 'k', 'linestyle', ':');
    annotation('textbox', 'position', [0.025 0.025 0.7 0.05], 'String', smFileList{f},...
        'linestyle', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
    drawnow;
    
    figure; imagesc(piWindow(1):1/1000:piWindow(2), 1:size(iscBetaPI,2),iscBetaPI', [-2 5]); colormap jet;
    annotation('textbox', 'position', [0.025 0.025 0.7 0.05], 'String', smFileList{f},...
        'linestyle', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
end

pkMtx = [([tdPI.PokeInIndex])', ([tdPI.PokeOutIndex])'];
figure; sp1 = subplot(3,1,1); plot(statMatrix(:,2)); sp2 = subplot(3,1,2); plot(statMatrix(:,betaCol)); sp3 = subplot(3,1,3); plot(zEnv); linkaxes([sp1 sp2 sp3], 'x'); hold on; plot(pkMtx', (ones(size(pkMtx))')*5, 'marker', '*', 'color', 'k');

%%
wlCorrVals = cell2mat(wlCorrVals');
tetIDs = tetIDs';
save('BetaWithdrawLatCorrSQR.mat', 'wlCorrVals', 'tetIDs');