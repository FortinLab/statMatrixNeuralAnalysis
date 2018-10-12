%% Spike-Phase-Power Plots
% Initial segment chopped from RMSpowerCoherenceEpocDetection.... maybe
% rework/clean up a bit.
%% Extract Info and bands and shit idk
origDir = cd;    
[fileDir] = uigetdir(origDir);
if fileDir==0
    disp('Analysis Cancelled')
    return
else  
    cd(fileDir)
end
dirContents = dir(fileDir);
fileNames = {dirContents.name};
tetFileLog = cellfun(@(a)~isempty(a), regexp(fileNames, '_T([0-9]*)_SM.mat'));
tetFiles = fileNames(tetFileLog)';
load(tetFiles{1}, 'statMatrixColIDs');
% Identify LFP VOLTAGE columns, identified by the fact that they do not
% contain the suffix '_HilbVals', i.e. their band name is their terminal
% item (determined by the $ in the regexp call below).
[bs, bf] = regexp(statMatrixColIDs, 'LFP_([A-Z,a-z]*)$');
% Identify the RAW column and remove it from the list since we only want
% the bandpass calls.
lfpRawLog = cellfun(@(a)~isempty(strfind(a, 'Raw')), statMatrixColIDs);
bs(lfpRawLog) = {[]};
bf(lfpRawLog) = {[]};
bandCols = find(cellfun(@(a)~isempty(a), bs));
bands = cell(length(bandCols),1);
% Now step through and extract the band names
for bnd = 1:length(bands)
    curCol = bandCols(bnd);
    bands{bnd} = statMatrixColIDs{curCol}(bs{curCol}+4:bf{curCol});
end
cd(fileDir);

%% Parameters
sampleRate = 1000;
windowSize = struct('Ripple', 0.008,...
    'HighGamma', 0.018,...
    'LowGamma', 0.026,...
    'Beta', 0.042,...
    'LowBeta', 0.1,...
    'Alpha', 0.127,... 
    'Theta', 0.22);

dirContents = dir(fileDir);
fileNames = {dirContents.name};
tetFileLog = cellfun(@(a)~isempty(a), regexp(fileNames, '_T([0-9]*)_SM.mat'));
tetFiles = fileNames(tetFileLog)';
saveYN = 1;

w = gausswin(21);
w = w/sum(w);
%% Analyzesing
for tet = 1:length(tetFiles)
    load(tetFiles{tet});
    uniColLog = cellfun(@(a)~isempty(a), regexp(statMatrixColIDs, '-U([0-9]*)'));
    numUnis = sum(uniColLog);
    if numUnis == 0
        continue
    else
        hilbTrace = cell(1,length(bands));
        rmsTrace = cell(1,length(bands));
        for bnd = 1:length(bands)
            windowNdxs = windowSize.(bands{bnd})*sampleRate;
            hlfWndo = windowNdxs/2;
            voltColLog = cellfun(@(a)~isempty(a), regexp(statMatrixColIDs,{[sprintf('_%s',bands{bnd}) '$']}));
            hilbColLog = cellfun(@(a)~isempty(a), regexp(statMatrixColIDs,{[sprintf('_%s',bands{bnd}) '_Hilb']}));
            hilbTrace{bnd} = statMatrix(:,hilbColLog);
            rmsTrace{bnd} = zscore(conv(sqrt(conv(statMatrix(:,voltColLog).^2, ones(windowNdxs,1)/windowNdxs, 'same')), w, 'same'));
        end
        uniSpots = find(uniColLog);
        for u = 1:numUnis
            curUniSpkLog = statMatrix(:,uniSpots(u))>0;
            figure('Name', statMatrixColIDs{uniSpots(u)});
            annotation('textbox', [0.05 0.9 0.9 0.1], 'String', sprintf('%s', statMatrixColIDs{uniSpots(u)}), 'linestyle', 'none', 'interpreter', 'none');
            for bnd = 1:length(bands)
                curUniPhasePref = hilbTrace{bnd}(curUniSpkLog);
                curUniPowerPref = rmsTrace{bnd}(curUniSpkLog);
                [sppPlot,Xedges,Yedges] = histcounts2(curUniPhasePref,curUniPowerPref,linspace(-pi,pi,13),0:0.5:10);
                subplot(1,length(bands),bnd)
                imagesc(-pi:0.1:(pi-0.1), 0:0.5:9.5, flipud(sppPlot'));
                set(gca, 'ytick', 0.5:1:9.5, 'yticklabel', 9:-1:0);
                title(bands{bnd});
                drawnow;                
            end
            if saveYN==1
                set(gcf, 'PaperOrientation', 'landscape');
                print('-fillpage', gcf, '-dpdf', sprintf('%s Spike-Phase-Power Plots', statMatrixColIDs{uniSpots(u)}));
            end
            close(gcf);
        end
    end
end
                
            
    
