freqs = 1:100;
for tet = 1:length(osmNeural)
    zData = zscore(osmNeural(tet).LFP);
%     zData = zscore(osmLFP(:,tet) - mean(osmLFP(:,2:end),2));
%     zData = osmLFP(:,tet) - mean(osmLFP(:,2:end),2);
    [s,fs,t] = spectrogram(zData, 1000, 900, freqs, 1000);
    ss = 10*log10(abs(s));
    zSS = ss;
    for f = 1:size(ss,1)
        zSS(f,:) = zscore(ss(f,:));
    end
    figure;
    sp1 = subplot(2,3,[1 2]);
    imagesc(t,fs,ss);
    hold on;
    scatter([osmBehavior.PortEntryTime], ones(1,length(osmBehavior))*20, '*k');
    scatter([osmBehavior.PortExitTime], ones(1,length(osmBehavior))*20, 'ok');
    set(gca, 'ydir', 'normal');
    title(osmNeural(tet).TetName);
    colorbar;
    
    sp2 = subplot(2,3,[4 5]);
    imagesc(t,fs,zSS, [-4 4]);
    hold on;
    scatter([osmBehavior.PortEntryTime], ones(1,length(osmBehavior))*20, '*k');
    scatter([osmBehavior.PortExitTime], ones(1,length(osmBehavior))*20, 'ok');
    set(gca, 'ydir', 'normal');
    colormap jet;
    linkaxes([sp1 sp2], 'xy');
    colorbar;
    
    subplot(2,3,3);
    periodogram(zData, [], freqs, 1000);
    
    subplot(2,3,6);
    plot(fs,median(zSS,2))
end
