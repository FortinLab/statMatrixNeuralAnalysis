%% Parameters
freqs = 1:100;
spectWin = 250;
spectOverlap = 200;
piWindow = [-0.5 0.5];
poWindow = [-0.5 0.5];

%%
spectNormMean = nan(length(freqs),length(osmNeural));
spectNormVar = nan(length(freqs),length(osmNeural));
for tet = 1:length(osmNeural)
    tempSsnSpect = spectrogram(osmNeural(tet).LFP, spectWin, spectOverlap, freqs, 1000);
    tempSsnSpect = 10*log10(abs(tempSsnSpect));
    spectNormMean(:,tet) = mean(tempSsnSpect,2);
    spectNormVar(:,tet) = std(tempSsnSpect,0,2);
end
for trl = 1:length(osmBehavior)
    fprintf('Trial #%i\n',trl);
    pokeInWindowLog = tsVect>osmBehavior(trl).PortEntryTime + piWindow(1) - (spectWin/2000) & tsVect<=osmBehavior(trl).PortEntryTime + piWindow(2) + (spectWin/2000);
    pokeOutWindowLog = tsVect>osmBehavior(trl).PortExitTime + poWindow(1) - (spectWin/2000) & tsVect<=osmBehavior(trl).PortExitTime + poWindow(2) + (spectWin/2000);    
    
    betaHilbPI = nan(sum(pokeInWindowLog), length(osmNeural));
    betaHilbPO = nan(sum(pokeOutWindowLog), length(osmNeural));
    osmBehavior(trl).ThetaPhaseBinPI = nan(length(0:spectWin-spectOverlap:diff(piWindow)*1000),length(osmNeural));
    osmBehavior(trl).ThetaPhaseBinPO = nan(length(0:spectWin-spectOverlap:diff(poWindow)*1000),length(osmNeural));

    
    osmBehavior(trl).PIspect = nan(length(freqs), length(0:spectWin-spectOverlap:diff(piWindow)*1000), length(osmNeural)); %#ok<*SAGROW>
    osmBehavior(trl).POspect = nan(length(freqs), length(0:spectWin-spectOverlap:diff(poWindow)*1000), length(osmNeural));
    for tet = 1:length(osmNeural)
        tetPokeInLFP = osmNeural(tet).LFP(pokeInWindowLog);
        tetPokeOutLFP = osmNeural(tet).LFP(pokeOutWindowLog);
        if tet == 1
            [tempPIspect, trlFreqs{trl}, ts] = spectrogram(tetPokeInLFP, spectWin, spectOverlap, freqs, 1000);
            tempTs = ts - (spectWin/2000);
            [tempPOspect, ~, ~] = spectrogram(tetPokeOutLFP, spectWin, spectOverlap, freqs, 1000);
        else
            [tempPIspect, ~, ~] = spectrogram(tetPokeInLFP, spectWin, spectOverlap, freqs, 1000);
            [tempPOspect, ~, ~] = spectrogram(tetPokeOutLFP, spectWin, spectOverlap, freqs, 1000);
        end
        osmBehavior(trl).PIspect(:,:,tet) = 10*log10(abs(tempPIspect));
        for t = 1:size(osmBehavior(trl).PIspect,2)
            osmBehavior(trl).PIspect(:,t,tet) = (osmBehavior(trl).PIspect(:,t,tet)-spectNormMean(:,tet))./spectNormVar(:,tet);
        end
        osmBehavior(trl).POspect(:,:,tet) = 10*log10(abs(tempPOspect));
        for t = 1:size(osmBehavior(trl).POspect,2)
            osmBehavior(trl).POspect(:,t,tet) = (osmBehavior(trl).POspect(:,t,tet)-spectNormMean(:,tet))./spectNormVar(:,tet);
        end
        [betaHilbPI(:,tet), ~] = PhaseFreqDetectAbbr(tetPokeInLFP, tsVect(pokeInWindowLog), 16, 32);
        [betaHilbPO(:,tet), ~] = PhaseFreqDetectAbbr(tetPokeOutLFP, tsVect(pokeOutWindowLog), 16, 32);
        [thetaHilbPI, ~] = PhaseFreqDetectAbbr(tetPokeInLFP, tsVect(pokeInWindowLog), 4, 12);
        piNndx = 0;
        for tBin = 0:spectWin-spectOverlap:diff(piWindow)*1000
            piNndx = piNndx+1;
            osmBehavior(trl).ThetaPhaseBinPI(piNndx,tet) = circ_mean(thetaHilbPI(tBin+1:tBin+(spectWin-spectOverlap)));
        end
        [thetaHilbPO, ~] = PhaseFreqDetectAbbr(tetPokeInLFP, tsVect(pokeOutWindowLog), 4, 12);
        poNndx = 0;
        for tBin = 0:spectWin-spectOverlap:diff(poWindow)*1000
            poNndx = poNndx+1;
            osmBehavior(trl).ThetaPhaseBinPO(poNndx,tet) = circ_mean(thetaHilbPO(tBin+1:tBin+(spectWin-spectOverlap)));
        end
    end   
    tempTrlCoherPI = nan(length(osmNeural), length(osmNeural), length(0:spectWin-spectOverlap:diff(piWindow)*1000));
    tempTrlCoherPO = nan(length(osmNeural), length(osmNeural), length(0:spectWin-spectOverlap:diff(piWindow)*1000));
    for tet1 = 1:length(osmNeural)
        for tet2 = 1:length(osmNeural)
            if tet1 == tet2
                continue;
            else
                tsNdx = 1;
                for ts = 1:(spectWin-spectOverlap):size(betaHilbPI,1)-(spectWin)+1
                    [tempTrlCoherPI(tet1, tet2, tsNdx), ~] = circ_corrcc(betaHilbPI(ts:ts+spectWin-1,tet1), betaHilbPI(ts:ts+spectWin-1,tet2));
                    [tempTrlCoherPO(tet1, tet2, tsNdx), ~] = circ_corrcc(betaHilbPO(ts:ts+spectWin-1,tet1), betaHilbPO(ts:ts+spectWin-1,tet2));
                    tsNdx = tsNdx+1;
                end
            end
        end        
    end    
    osmBehavior(trl).PIcoher = tempTrlCoherPI;
    osmBehavior(trl).PIcoherMean = nanmean(tempTrlCoherPI,3);
    osmBehavior(trl).POcoher = tempTrlCoherPO;
    osmBehavior(trl).POcoherMean = nanmean(tempTrlCoherPO,3);
end

%%
isTrlLog = ([osmBehavior.TrialPosition] - [osmBehavior.TrialOdor])==0;
perfTrlLog = [osmBehavior.PerformanceLog]==1;
iscTrials = osmBehavior(isTrlLog & perfTrlLog);
oscTrials = osmBehavior(~isTrlLog & perfTrlLog);

for tet = 1:length(osmNeural)
    iscPIspect = median(cell2mat(reshape(cellfun(@(a)a(:,:,tet), {iscTrials.PIspect}, 'uniformoutput', 0), [1,1,length(iscTrials)])),3);
    oscPIspect = median(cell2mat(reshape(cellfun(@(a)a(:,:,tet), {oscTrials.PIspect}, 'uniformoutput', 0), [1,1,length(oscTrials)])),3);
    iscPOspect = median(cell2mat(reshape(cellfun(@(a)a(:,:,tet), {iscTrials.POspect}, 'uniformoutput', 0), [1,1,length(iscTrials)])),3);
    oscPOspect = median(cell2mat(reshape(cellfun(@(a)a(:,:,tet), {oscTrials.POspect}, 'uniformoutput', 0), [1,1,length(oscTrials)])),3);
    figure;
    annotation('textbox', 'Position',[0.05 0.95 0.9 0.05],...
        'String', osmNeural(tet).TetName, 'HorizontalAlignment', 'Left',...
        'Fontweight', 'bold',...
        'VerticalAlignment', 'bottom', 'FitHeightToText','off',...
        'LineStyle','none', 'Interpreter', 'none');
    subplot(2,2,1);
    imagesc((0:spectWin-spectOverlap:diff(piWindow)*1000)-500, freqs, iscPIspect, [-2 2]);
    title('ISC Poke In Aligned');
    set(gca, 'ydir', 'normal');
    subplot(2,2,2);
    imagesc((0:spectWin-spectOverlap:diff(poWindow)*1000)-500, freqs, iscPOspect, [-2 2]);
    title('ISC Poke Out Aligned');
    set(gca, 'ydir', 'normal');
    subplot(2,2,3);
    imagesc((0:spectWin-spectOverlap:diff(piWindow)*1000)-500, freqs, oscPIspect, [-2 2]);
    title('ISC Poke In Aligned');
    set(gca, 'ydir', 'normal');
    subplot(2,2,4);
    imagesc((0:spectWin-spectOverlap:diff(poWindow)*1000)-500, freqs, oscPOspect, [-2 2]);
    title('ISC Poke Out Aligned');
    set(gca, 'ydir', 'normal');
end

figure;
iscPIcoher = median(cell2mat(reshape({iscTrials.PIcoherMean}, [1,1,length(iscTrials)])),3);
oscPIcoher = median(cell2mat(reshape({oscTrials.PIcoherMean}, [1,1,length(oscTrials)])),3);
iscPOcoher = median(cell2mat(reshape({iscTrials.POcoherMean}, [1,1,length(iscTrials)])),3);
oscPOcoher = median(cell2mat(reshape({oscTrials.POcoherMean}, [1,1,length(oscTrials)])),3);
subplot(2,2,1);
imagesc(1:length(osmNeural), 1:length(osmNeural), iscPIcoher, [-1 1]);
set(gca, 'xtick', 1:length(osmNeural), 'xticklabel', {osmNeural.TetName}, 'ytick', 1:length(osmNeural), 'yticklabel', {osmNeural.TetName});
title('ISC Beta Coherence (PI)');

subplot(2,2,2);
imagesc(1:length(osmNeural), 1:length(osmNeural), iscPOcoher, [-1 1]);
set(gca, 'xtick', 1:length(osmNeural), 'xticklabel', {osmNeural.TetName}, 'ytick', 1:length(osmNeural), 'yticklabel', {osmNeural.TetName});
title('ISC Beta Coherence (PO)');

subplot(2,2,3);
imagesc(1:length(osmNeural), 1:length(osmNeural), oscPIcoher, [-1 1]);
set(gca, 'xtick', 1:length(osmNeural), 'xticklabel', {osmNeural.TetName}, 'ytick', 1:length(osmNeural), 'yticklabel', {osmNeural.TetName});
title('OSC Beta Coherence (PI)');

subplot(2,2,4);
imagesc(1:length(osmNeural), 1:length(osmNeural), oscPOcoher, [-1 1]);
set(gca, 'xtick', 1:length(osmNeural), 'xticklabel', {osmNeural.TetName}, 'ytick', 1:length(osmNeural), 'yticklabel', {osmNeural.TetName});
title('OSC Beta Coherence (PI)');
colormap jet

%%
phaseBins = -pi:pi/5:pi;
for tet = 1:length(osmNeural)
    tpbPI = cell2mat(cellfun(@(a)a(:,tet),{iscTrials.ThetaPhaseBinPI}, 'uniformoutput', 0));
    tpbPO = cell2mat(cellfun(@(a)a(:,tet),{iscTrials.ThetaPhaseBinPO}, 'uniformoutput', 0));
    betaPowerPI = cell2mat(cellfun(@(a)mean(a(16:32,:,tet),1)', {iscTrials.PIspect}, 'uniformoutput', 0));
    betaPowerPO = cell2mat(cellfun(@(a)mean(a(16:32,:,tet),1)', {iscTrials.POspect}, 'uniformoutput', 0));
    
    piThetaBetaMod = nan(length(phaseBins)-1, length(0:spectWin-spectOverlap:diff(piWindow)*1000));
    for t = 1:length(0:spectWin-spectOverlap:diff(piWindow)*1000)
        curTPB = tpbPI(t,:);
        for p = 1:length(phaseBins)-1
            tpbLog = curTPB>=phaseBins(p) & curTPB<phaseBins(p+1);
            piThetaBetaMod(p,t) = mean(betaPowerPI(tpbLog));
        end
    end
            
    poThetaBetaMod = nan(length(phaseBins)-1, length(0:spectWin-spectOverlap:diff(poWindow)*1000));
    for t = 1:length(0:spectWin-spectOverlap:diff(piWindow)*1000)
        curTPB = tpbPO(t,:);
        for p = 1:length(phaseBins)-1
            tpbLog = curTPB>=phaseBins(p) & curTPB<phaseBins(p+1);
            poThetaBetaMod(p,t) = mean(betaPowerPO(tpbLog));
        end
    end
    figure;
    annotation('textbox', 'Position',[0.05 0.95 0.9 0.05],...
        'String', osmNeural(tet).TetName, 'HorizontalAlignment', 'Left',...
        'Fontweight', 'bold',...
        'VerticalAlignment', 'bottom', 'FitHeightToText','off',...
        'LineStyle','none', 'Interpreter', 'none');
    subplot(1,2,1)
    imagesc((0:spectWin-spectOverlap:diff(piWindow)*1000)-500, phaseBins(1:end-1)+mode(diff(phaseBins)/2), piThetaBetaMod, [-1 1]);
    title('Theta-Beta Phase Amplitude Coupling (Poke In)');
    subplot(1,2,2)
    imagesc((0:spectWin-spectOverlap:diff(poWindow)*1000)-500, phaseBins(1:end-1)+mode(diff(phaseBins)/2), poThetaBetaMod, [-1 1]);
    title('Theta-Beta Phase Amplitude Coupling (Poke In)');
end
    

                