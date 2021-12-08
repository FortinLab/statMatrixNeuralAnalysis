function PlotRipFeatsByOdor(rips, trialRips)

%% Extract Logicals
perfLog = rips.TrialInfo.Perf;
isLog = rips.TrialInfo.TransDist==0;

%% Plot Data
figure;
% Plot Odor A
aLog = perfLog & isLog & rips.TrialInfo.OdorVect==1;
aP = subplot(4,5,1);
    RipLikelyPlot(trialRips.Events, aLog);
    ylabel(aP, '% Trials w/Ripple');
    title(aP, 'Odor A');
aD = subplot(4,5,6);
    RipPlotVar(trialRips.Duration, aLog);
    ylabel(aD, 'Duration');
aS = subplot(4,5,11);
    RipPlotVar(trialRips.Synchrony, aLog);
    ylabel(aS, 'Synchrony');
aE = subplot(4,5,16);
    RipPlotVar(trialRips.EnsembleAct, aLog);
    ylabel(aE, '% Units Active');
    set(aE, 'xticklabels', [{'Pre-Trial'}, {'Trial'}, {'Post-Trial'}], 'box', 'off');
    
% Plot Odor B
bLog = perfLog & isLog & rips.TrialInfo.OdorVect==2;
bP = subplot(4,5,2);
    RipLikelyPlot(trialRips.Events, bLog);
    ylabel(bP, '% Trials w/Ripple');
    title(bP, 'Odor B');
bD = subplot(4,5,7);
    RipPlotVar(trialRips.Duration, bLog);
    ylabel(bD, 'Duration');
bS = subplot(4,5,12);
    RipPlotVar(trialRips.Synchrony, bLog);
    ylabel(bS, 'Synchrony');
bE = subplot(4,5,17);
    RipPlotVar(trialRips.EnsembleAct, bLog);
    ylabel(bE, '% Units Active');
    set(bE, 'xticklabels', [{'Pre-Trial'}, {'Trial'}, {'Post-Trial'}], 'box', 'off');
    
% Plot Odor C
cLog = perfLog & isLog & rips.TrialInfo.OdorVect==3;
cP = subplot(4,5,3);
    RipLikelyPlot(trialRips.Events, cLog);
    ylabel(cP, '% Trials w/Ripple');
    title(cP, 'Odor C');
cD = subplot(4,5,8);
    RipPlotVar(trialRips.Duration, cLog);
    ylabel(cD, 'Duration');
cS = subplot(4,5,13);
    RipPlotVar(trialRips.Synchrony, cLog);
    ylabel(cS, 'Synchrony');
cE = subplot(4,5,18);
    RipPlotVar(trialRips.EnsembleAct, cLog);
    ylabel(cE, '% Units Active');
    set(cE, 'xticklabels', [{'Pre-Trial'}, {'Trial'}, {'Post-Trial'}], 'box', 'off');
    
% Plot Odor D
dLog = perfLog & isLog & rips.TrialInfo.OdorVect==4;
dP = subplot(4,5,4);
    RipLikelyPlot(trialRips.Events, dLog);
    ylabel(dP, '% Trials w/Ripple');
    title(dP, 'Odor D');
dD = subplot(4,5,9);
    RipPlotVar(trialRips.Duration, dLog);
    ylabel(dD, 'Duration');
dS = subplot(4,5,14);
    RipPlotVar(trialRips.Synchrony, dLog);
    ylabel(dS, 'Synchrony');
dE = subplot(4,5,19);
    RipPlotVar(trialRips.EnsembleAct, dLog);
    ylabel(dE, '% Units Active');
    set(dE, 'xticklabels', [{'Pre-Trial'}, {'Trial'}, {'Post-Trial'}], 'box', 'off');
    
% Plot Odor E
eLog = perfLog & isLog & rips.TrialInfo.OdorVect==5;
eP = subplot(4,5,5);
    RipLikelyPlot(trialRips.Events, eLog);
    ylabel(eP, '% Trials w/Ripple');
    title(eP, 'Odor E');
eD = subplot(4,5,10);
    RipPlotVar(trialRips.Duration, eLog);
    ylabel(eD, 'Duration');
eS = subplot(4,5,15);
    RipPlotVar(trialRips.Synchrony, eLog);
    ylabel(eS, 'Synchrony');
eE = subplot(4,5,20);
    RipPlotVar(trialRips.EnsembleAct, eLog);
    ylabel(eE, '% Units Active');
    set(eE, 'xticklabels', [{'Pre-Trial'}, {'Trial'}, {'Post-Trial'}], 'box', 'off');

linkaxes([aP bP cP dP eP], 'xy');
linkaxes([aD bD cD dD eD], 'xy');
linkaxes([aS bS cS dS eS], 'xy');
linkaxes([aE bE cE dE eE], 'xy');
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
