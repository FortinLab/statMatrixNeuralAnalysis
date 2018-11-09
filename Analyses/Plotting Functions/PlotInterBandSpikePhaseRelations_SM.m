function PlotInterBandSpikePhaseRelations_SM(statMatrix, statMatrixColIDs, saveYN)
%% Identify Input Variables
% Identify the LFP bands filtered
phaseValsLog = cellfun(@(a)~isempty(a), regexp(statMatrixColIDs, '_HilbVals\>'));
rawChanLog = cellfun(@(a)~isempty(a), strfind(statMatrixColIDs, 'Raw'));
phaseValCols = find(phaseValsLog & ~rawChanLog);
lfpBands = cellfun(@(a)a{3}, cellfun(@(a)strsplit(a, '_'), statMatrixColIDs(phaseValCols), 'uniformoutput', 0), 'uniformoutput', 0);

% Identify units
uniLog = cellfun(@(a)~isempty(a), strfind(statMatrixColIDs, '-U'));
unitIDs = statMatrixColIDs(uniLog);
unitCols = find(uniLog);

%% Define Plotting Constants
phaseBins = -pi:pi/8:pi;
%%
for uni = 1:sum(uniLog)
    figure
    colPlots = nan(length(lfpBands));
    cLimVals = cell(length(lfpBands));
    annotation('textbox', [0.05 0.9 0.9 0.1], 'String', sprintf('%s Session Spiking Inter-Band Phase Relations', unitIDs{uni}), 'linestyle', 'none');
    curUniSpkLog = statMatrix(:,unitCols(uni))>=1;
    for bnd1 = 1:length(lfpBands)
        bnd1Phase = statMatrix(curUniSpkLog,phaseValCols(bnd1));
        for bnd2 = 1:length(lfpBands)
            bnd2Phase = statMatrix(curUniSpkLog,phaseValCols(bnd2));
            sp = sub2ind([length(lfpBands),length(lfpBands)], bnd1, bnd2);
            colPlots(sp) = subplot(length(lfpBands), length(lfpBands), sp); 
            contourf(phaseBins(1:end-1),phaseBins(1:end-1),histcounts2(bnd1Phase, bnd2Phase, phaseBins, phaseBins)', 20, 'linestyle', 'none');
            if bnd1~=bnd2
                cLimVals{sp} = get(colPlots(sp), 'CLim');
            end
            set(colPlots(sp), 'xticklabel', [], 'yticklabel', [], 'DataAspectRatio',[1 1 1],'Layer','top');
            axis on;            
            xlabel(sprintf('%s Phase', lfpBands{bnd1}));
            ylabel(sprintf('%s Phase', lfpBands{bnd2}));
        end
    end
    commonClims = [min(min(cell2mat(cLimVals(:)))), max(max(cell2mat(cLimVals(:))))*.75];
    for sp = 1:length(lfpBands)^2
        set(colPlots(sp), 'clim', commonClims);
    end
    colormap jet
    if saveYN==1
        set(gcf, 'PaperOrientation', 'landscape');
        print('-fillpage', gcf, '-dpdf', sprintf('%s Spiking Inter-Band Phase Relations', unitIDs{uni}));
    end
    drawnow
end
