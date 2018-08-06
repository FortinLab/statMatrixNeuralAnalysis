z = horzcat(statDiffPerTetPERF{:})
xTicks = (eventWindow(1):pehBinSize:(eventWindow(2)-pehBinSize))+(pehBinSize/2);

pokeIn = [z.PokeIn];
pokeOut = [z.PokeOut];

pokeInAll = cell2mat({pokeIn.All}');
pokeOutAll = cell2mat({pokeOut.All}');
pokeInAllMean = nanmean(pokeInAll);
pokeOutAllMean = nanmean(pokeOutAll);
pokeInAllSD = nanstd(pokeInAll,0,1)./(sum(~isnan(pokeInAll))-1);
pokeOutAllSD = nanstd(pokeOutAll,0,1)./(sum(~isnan(pokeOutAll))-1);

pokeInCorr = cell2mat({pokeIn.Correct}');
pokeOutCorr = cell2mat({pokeOut.Correct}');
pokeInCorrMean = nanmean(pokeInCorr);
pokeOutCorrMean = nanmean(pokeOutCorr);
pokeInCorrSD = nanstd(pokeInCorr,0,1)./(sum(~isnan(pokeInCorr))-1);
pokeOutCorrSD = nanstd(pokeOutCorr,0,1)./(sum(~isnan(pokeOutCorr))-1);

pokeInInCorr = cell2mat({pokeIn.Incorrect}');
pokeOutInCorr = cell2mat({pokeOut.Incorrect}');
pokeInInCorrMean = nanmean(pokeInInCorr);
pokeOutInCorrMean = nanmean(pokeOutInCorr);
pokeInInCorrSD = nanstd(pokeInInCorr,0,1)./(sum(~isnan(pokeInInCorr))-1);
pokeOutInCorrSD = nanstd(pokeOutInCorr,0,1)./(sum(~isnan(pokeOutInCorr))-1);

figure; 
ax1 = subplot(2,1,1);
piAll = PlotLineAndFilledError(xTicks, pokeInAllMean, pokeInAllSD, 'black');
hold on
piCorr = PlotLineAndFilledError(xTicks, pokeInCorrMean, pokeInCorrSD, 'red');
piInCorr = PlotLineAndFilledError(xTicks, pokeInInCorrMean, pokeInInCorrSD, 'Blue');
legend([piAll piCorr piInCorr], 'All Trials', 'Correct Trials', 'Incorrect Trials');
line(get(ax1, 'xlim'), [0 0], 'linestyle', ':', 'linewidth', 1.5, 'color','black');
axis tight
title('Poke In Aligned F-Ratio Differences');
xlabel('Time Relative to Poke Initiation');
ylabel({'F-Ratio Difference', '(Position - Odor)'});

ax2 = subplot(2,1,2);
poAll = PlotLineAndFilledError(xTicks, pokeOutAllMean, pokeOutAllSD, 'black');
hold on
poCorr = PlotLineAndFilledError(xTicks, pokeOutCorrMean, pokeOutCorrSD, 'red');
poInCorr = PlotLineAndFilledError(xTicks, pokeOutInCorrMean, pokeOutInCorrSD, 'Blue');
line(get(ax1, 'xlim'), [0 0], 'linestyle', ':', 'linewidth', 1.5, 'color','black');
axis tight
title('Poke Out Aligned F-Ratio Differences');
xlabel('Time Relative to Poke Withdrawal');
ylabel({'F-Ratio Difference', '(Position - Odor)'});

linkaxes([ax1 ax2], 'xy')
