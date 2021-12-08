function PlotRipFeatsByTrlType(rips, trialRips)

%% Extract Logicals
perfLog = rips.TrialInfo.Perf;
isLog = rips.TrialInfo.TransDist==0;

%% Plot Data
% Plot InSeq Correct
iscLog = perfLog & isLog;
figure;
iscP = subplot(4,4,1);
    RipLikelyPlot(trialRips.Events, iscLog);
    ylabel(iscP, '% Trials w/Ripple');
    title(iscP, 'InSeq Correct Trials');
iscD = subplot(4,4,5);
    RipPlotVar(trialRips.Duration, iscLog);
    ylabel(iscD, 'Duration');
iscS = subplot(4,4,9);
    RipPlotVar(trialRips.Synchrony, iscLog);
    ylabel(iscS, 'Synchrony');
iscE = subplot(4,4,13);
    RipPlotVar(trialRips.EnsembleAct, iscLog);
    ylabel(iscE, '% Units Active');
    set(gca, 'xticklabels', [{'Pre-Trial'}, {'Trial'}, {'Post-Trial'}], 'box', 'off');
    
% Plot InSeq Incorrect
isiLog = ~perfLog & isLog;
isiP = subplot(4,4,2);
    RipLikelyPlot(trialRips.Events, isiLog);
    ylabel(isiP, '% Trials w/Ripple');
    title(isiP, 'InSeq Incorrect Trials');
isiD = subplot(4,4,6);
    RipPlotVar(trialRips.Duration, isiLog);
    ylabel(isiD, 'Duration');
isiS = subplot(4,4,10);
    RipPlotVar(trialRips.Synchrony, isiLog);
    ylabel(isiS, 'Synchrony');
isiE = subplot(4,4,14);
    RipPlotVar(trialRips.EnsembleAct, isiLog);
    ylabel(isiE, '% Units Active');
    set(gca, 'xticklabels', [{'Pre-Trial'}, {'Trial'}, {'Post-Trial'}], 'box', 'off');

% Plot OutSeq Incorrect
osiLog = ~perfLog & ~isLog;
osiP = subplot(4,4,3);
    RipLikelyPlot(trialRips.Events, osiLog);
    title(osiP, 'OutSeq Incorrect Trials');
osiD = subplot(4,4,7);
    RipPlotVar(trialRips.Duration, osiLog);
osiS = subplot(4,4,11);
    RipPlotVar(trialRips.Synchrony, osiLog);
osiE = subplot(4,4,15);
    RipPlotVar(trialRips.EnsembleAct, osiLog);
    set(gca, 'xticklabels', [{'Pre-Trial'}, {'Trial'}, {'Post-Trial'}], 'box', 'off');

% Plot OutSeq Correct
oscLog = perfLog & ~isLog;
oscP = subplot(4,4,4);
    RipLikelyPlot(trialRips.Events, oscLog);
    title(oscP, 'OutSeq Correct Trials');
oscD = subplot(4,4,8);
    RipPlotVar(trialRips.Duration, oscLog);
oscS = subplot(4,4,12);
    RipPlotVar(trialRips.Synchrony, oscLog);
oscE = subplot(4,4,16);
    RipPlotVar(trialRips.EnsembleAct, oscLog);
   set(gca, 'xticklabels', [{'Pre-Trial'}, {'Trial'}, {'Post-Trial'}], 'box', 'off');


linkaxes([iscP isiP oscP osiP], 'xy');
linkaxes([iscD isiD oscD osiD], 'xy');
linkaxes([iscS isiS oscS osiS], 'xy');
linkaxes([iscE isiE oscE osiE], 'xy');
end


%% Plotting Functions
function RipLikelyPlot(trialRips, log)
    prcntTrlsWrips = mean(cellfun(@(a)~isempty(a),trialRips(log,:)));
    bar(prcntTrlsWrips);
    set(gca, 'xticklabels', [], 'box', 'off');
end

function RipPlotVar(trialRipsParam, log)
    BarPlotErrorbars([mean(cell2mat(trialRipsParam(log,1))),mean(cell2mat(trialRipsParam(log,2))),mean(cell2mat(trialRipsParam(log,3)))],...
        [std(cell2mat(trialRipsParam(log,1))),std(cell2mat(trialRipsParam(log,2))),std(cell2mat(trialRipsParam(log,3)))],[],[]);
    set(gca, 'xticklabels', [], 'box', 'off');
end