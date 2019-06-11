[unitEpoch, unitIDs, lfpEpoch, lfpIDs, trialTimeBins, eventTimeBins, trialInfo] = EpochExtraction_SM('PokeIn', -0.5, 1.5, 'lfpBand', 'Theta', 'lfpData', 'Both');

% Parse trials and do selection etc
iscLog = trialInfo(:,1)==1 & trialInfo(:,2)==1;

% Extract the reference LFP values
lfpRef = rad2deg(lfpEpoch(:,:,strcmp('T14_LFP_Theta_HilbVals', lfpIDs)));
timeMtx = repmat(eventTimeBins', [size(lfpRef,1), 1]);
for u = 1:size(unitEpoch,3)
    curUni = logical(unitEpoch(:,:,u));
    tempUni = curUni(iscLog,:);
    tempLFP = lfpRef(iscLog,:);
    tempTime = timeMtx(iscLog,:);
    figure;
    subplot(6,1,1)
    scatter(tempTime(tempUni), tempLFP(tempUni), '.k');
    set(gca, 'ylim', [-180 180], 'xlim', [-0.5 1.5]);
    title('All InSeq Correct Trials');
    
    for p = 1:5
        posLog = trialInfo(:,3)==p;
        tempUni = curUni(iscLog & posLog,:);
        tempLFP = lfpRef(iscLog & posLog,:);
        tempTime = timeMtx(iscLog & posLog,:);
        subplot(6,1,p+1)
        scatter(tempTime(tempUni), tempLFP(tempUni), '.k');
        set(gca, 'ylim', [-180 180], 'xlim', [-0.5 1.5]);
        title(sprintf('Position %i InSeq Correct', p));
    end
    annotation('textbox', [0.05 0.9 0.9 0.1], 'String', sprintf('%s', unitIDs{u}), 'linestyle', 'none', 'interpreter', 'none');
    set(gcf, 'PaperOrientation', 'landscape');
    print('-fillpage', gcf, '-dpdf', [unitIDs{u} '_Spike_Phase_Scatter.pdf']);
    close gcf
end
    