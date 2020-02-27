function PlotNearTrialRipStats(tRips)
trialRips = tRips.Events;
trialRipsDur = tRips.Duration;
trialRipsSync = tRips.Synchrony;
trialRipsNsmbl = tRips.EnsembleAct;

figure;
prcntTrlsWrips = mean(cellfun(@(a)~isempty(a),trialRips));
subplot(4,1,1);
bar(prcntTrlsWrips);
set(gca, 'xticklabels', [], 'box', 'off');
ylabel('% Trials w/Ripple');
subplot(4,1,2);
BarPlotErrorbars([mean(cell2mat(trialRipsDur(:,1))),mean(cell2mat(trialRipsDur(:,2))),mean(cell2mat(trialRipsDur(:,3)))],...
    [std(cell2mat(trialRipsDur(:,1))),std(cell2mat(trialRipsDur(:,2))),std(cell2mat(trialRipsDur(:,3)))],[],[]);
set(gca, 'xticklabels', [], 'box', 'off');
ylabel('Duration (ms)');
subplot(4,1,3);
BarPlotErrorbars([mean(cell2mat(trialRipsSync(:,1))),mean(cell2mat(trialRipsSync(:,2))),mean(cell2mat(trialRipsSync(:,3)))],...
    [std(cell2mat(trialRipsSync(:,1))),std(cell2mat(trialRipsSync(:,2))),std(cell2mat(trialRipsSync(:,3)))],[],[]);
set(gca, 'xticklabels', [{'Pre-Trial'}, {'Trial'}, {'Post-Trial'}], 'box', 'off');
ylabel('Synchrony');
subplot(4,1,4);
BarPlotErrorbars([mean(cell2mat(trialRipsNsmbl(:,1))),mean(cell2mat(trialRipsNsmbl(:,2))),mean(cell2mat(trialRipsNsmbl(:,3)))],...
    [std(cell2mat(trialRipsNsmbl(:,1))),std(cell2mat(trialRipsNsmbl(:,2))),std(cell2mat(trialRipsNsmbl(:,3)))],[],[]);
set(gca, 'xticklabels', [{'Pre-Trial'}, {'Trial'}, {'Post-Trial'}], 'box', 'off');
ylabel('% Active Units');
end