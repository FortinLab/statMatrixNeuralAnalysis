function PlotRipCountsByEvent(rips)
%%
pokeInNdxs = rips.TrialInfo.TrialPokes(:,1);
pokeOutNdxs = rips.TrialInfo.TrialPokes(:,2);
rwdNdxs = rips.TrialInfo.TrialRewards;

if length(pokeInNdxs)~=length(rwdNdxs)
    rwdNdxs(end:end+(length(pokeInNdxs)-length(rwdNdxs)))=nan;
end
%%
piDiffs = cell(length(pokeInNdxs),1);
poDiffs = piDiffs;
rdDiffs = piDiffs;
for trl = 1:size(rips.TrialInfo.TrialPokes,1)
    curPIdiff = rips.Ripples.Events(:,1) - pokeInNdxs(trl);
    piDiffs{trl} = curPIdiff(curPIdiff>=-500 & curPIdiff<=500);
    
    curPOdiff = rips.Ripples.Events(:,1) - pokeOutNdxs(trl);
    poDiffs{trl} = curPOdiff(curPOdiff>=-500 & curPOdiff<=500);
    
    curRDdiff = rips.Ripples.Events(:,1) - rwdNdxs(trl);
    rdDiffs{trl} = curRDdiff(curRDdiff>=-500 & curRDdiff<=500);
end  
    
%%
figure;
subplot(3,1,1)
histogram(cell2mat(piDiffs), -500:10:500);
title('Ripples Relative to Poke In');
subplot(3,1,2)
histogram(cell2mat(poDiffs), -500:10:500);
title('Ripples Relative to Poke Out');
subplot(3,1,3)
histogram(cell2mat(rdDiffs), -500:10:500);
title('Ripples Relative to Reward Delivery');

