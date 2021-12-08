function PlotRipFeatCorr(param1, param2, param1Title, param2Title)
%% Plot Duration vs Synchrony Correlation
figure;
p1Hist = subplot(4,4,[1,5,9]);
histogram(param1, 30, 'Orientation', 'Horizontal');
set(gca, 'XDir', 'reverse');
line(get(gca, 'xlim'), repmat(mean(param1), [1,2]), 'color', 'r', 'linewidth', 2);
line(get(gca, 'xlim'), repmat(median(param1), [1,2]),'color', 'r', 'linestyle', '--', 'linewidth', 2);
if mean(param1)>median(param1)
    set(gca, 'ytick', [median(param1), mean(param1)], 'yticklabel', [{'Median'}, {'Mean'}]);
else
    set(gca, 'ytick', [mean(param1), median(param1)], 'yticklabel', [{'Mean'}, {'Median'}]);
end
xlabel('Count');
ylabel([{['\bf' param1Title '\rm']}; sprintf('Mean: %0.2f; Median: %0.2f', mean(param1), median(param1))]);

corre = subplot(4,4,[2:4,6:8,10:12]);
corrScatPlot(param2, param1, [], [], []);
title(sprintf('Session Wide Ripple %s vs %s', param1Title, param2Title));
set(gca, 'xticklabel', [], 'yticklabel', []);

p2Hist = subplot(4,4,14:16);
histogram(param2, 30);
title([{['\bf' param2Title '\rm']}; sprintf('Mean: %0.2f; Median: %0.2f', mean(param2), median(param2))]);
ylabel('Count');

linkaxes([p1Hist, corre], 'y');
linkaxes([p2Hist, corre], 'x');