odor = [bmts.Odor];
pos = [bmts.Position];
performance = [bmts.Performance];
isLog = [bmts.TranspositionDistance]==0;

for t = 1:length(ensembleMatrixColIDs)-1
    figure;
    sp = nan(1,4);
    for o = 1:4
        curOdorFam = trialOrgData(odor==o & isLog==1 & performance==1);
        curFam = mean(cell2mat(cellfun(@(a)conv(a(:,t),gausswin(200)./max(gausswin(200)),'same'),curOdorFam, 'uniformoutput',0)),2);
        curOdorNov = trialOrgData(odor==o+10 & isLog==1 & performance==1);
        curNov = mean(cell2mat(cellfun(@(a)conv(a(:,t),gausswin(200)./max(gausswin(200)),'same'),curOdorNov, 'uniformoutput',0)),2);
        sp(o) = subplot(1,4,o);
        plot(-500:1750, curFam, '-k');
        hold on;
        plot(-500:1750, curNov, '--k');
        title(sprintf('Position %i', o));
        xlabel('Time From Poke In (ms)')
    end
    linkaxes(sp, 'xy');
    legend('Familiar', 'Novel')
    for o = 1:4
        plot(sp(o), [0 0], get(sp(o), 'ylim'), '-k');
    end
    axis tight;
    annotation(gcf,'textbox', [0 0.93 1 0.05],'String', ensembleMatrixColIDs(t+1),...
    'FontWeight','bold', 'FontSize',12, 'edgecolor', 'none', 'horizontalalignment', 'left');
end
        
    
