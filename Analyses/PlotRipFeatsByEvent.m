function PlotRipFeatsByEvent(rips, param)
%%
pokeInNdxs = rips.TrialInfo.TrialPokes(:,1);
pokeOutNdxs = rips.TrialInfo.TrialPokes(:,2);
rwdNdxs = rips.TrialInfo.TrialRewards;

if length(pokeInNdxs)~=length(rwdNdxs)
    rwdNdxs(end:end+(length(pokeInNdxs)-length(rwdNdxs)))=nan;
end

%%
switch param
    case 'Duration'
        feat = rips.Ripples.Duration;
    case 'Synchrony'
        feat = rips.Ripples.Synchrony;
    case 'Spiking'
        feat = rips.Ripples.EnsembleActivity;
    case 'Power'
        feat = rips.Ripples.MaxPower;
    case 'MaxFreq'
        feat = rips.Ripples.MaxPowerFrequency;
end
        
%%
piDiffs = cell(length(pokeInNdxs),1);
piFeat = piDiffs;
poDiffs = piDiffs;
poFeat = piDiffs;
rdDiffs = piDiffs;
rdFeat = piDiffs;
for trl = 1:size(rips.TrialInfo.TrialPokes,1)
    curPIdiff = rips.Ripples.Events(:,1) - pokeInNdxs(trl);
    piLog = curPIdiff>=-500 & curPIdiff<=500;
    piDiffs{trl} = curPIdiff(piLog);
    piFeat{trl} = feat(piLog);
    
    curPOdiff = rips.Ripples.Events(:,1) - pokeOutNdxs(trl);
    poLog = curPOdiff>=-500 & curPOdiff<=500;
    poDiffs{trl} = curPOdiff(poLog);
    poFeat{trl} = feat(poLog);
    
    curRDdiff = rips.Ripples.Events(:,1) - rwdNdxs(trl);
    rdLog = curRDdiff>=-500 & curRDdiff<=500;
    rdDiffs{trl} = curRDdiff(rdLog);
    rdFeat{trl} = feat(rdLog);
end  
    
%%
figure;
sp1 = subplot(3,1,1);
scatter(cell2mat(piDiffs), cell2mat(piFeat), '*k');
title(sprintf('Ripple %s Relative to Poke In', param));
sp2 = subplot(3,1,2);
scatter(cell2mat(poDiffs), cell2mat(poFeat), '*k');
title(sprintf('Ripple %s Relative to Poke Out', param));
sp3 = subplot(3,1,3);
scatter(cell2mat(rdDiffs), cell2mat(rdFeat), '*k');
title(sprintf('Ripple %s Relative to Reward Delivery', param));
linkaxes([sp1 sp2 sp3], 'xy');
set(sp3, 'xlim', [-500 500]);