
% Phase-Power Coupling
%
%

%% Load Data
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
lfpRawLog = cellfun(@(a)~isempty(strfind(a, 'Raw')), statMatrixColIDs); %#ok<*STREMP>
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

%% Determine LFP band Parameters
sampleRate = 1000;
if sum(strcmp(bands, 'Alpha'))==0
%     windowSize = struct('Ripple', 0.008,...
%         'HighGamma', 0.018,...
%         'LowGamma', 0.026,...
%         'Beta', 0.042,...
%         'LowBeta', 0.1,...
%         'Theta', 0.22);
    windowSize = struct('HighGamma', 0.018,...
        'LowGamma', 0.026,...
        'Beta', 0.042,...
        'LowBeta', 0.1,...
        'Theta', 0.22);
else
%     windowSize = struct('Ripple', 0.008,...
%         'HighGamma', 0.018,...
%         'LowGamma', 0.026,...
%         'Beta', 0.042,...
%         'LowBeta', 0.1,...
%         'Alpha', 0.127,...
%         'Theta', 0.22);
    windowSize = struct('HighGamma', 0.018,...
        'LowGamma', 0.026,...
        'Beta', 0.042,...
        'LowBeta', 0.1,...
        'Alpha', 0.127,...
        'Theta', 0.22);
end

w = gausswin(21);
w = w/sum(w);
%% Identify/Extract LFP data
tetHilbVals = cell(1,length(tetFiles));
tetRmsVals = cell(1,length(tetFiles));
for fl = 1:length(tetFiles)
    load(tetFiles{fl});
    hilbVals = cell(1,length(bands));
    rmsVals = cell(1,length(bands));
    for b = 1:length(bands)
        if sum(strcmp(fieldnames(windowSize), bands{b}))>=1
            windowNdxs = windowSize.(bands{b})*sampleRate;
            voltColLog = cellfun(@(a)~isempty(a), regexp(statMatrixColIDs,{[sprintf('_%s',bands{b}) '$']}));
            hilbColLog = cellfun(@(a)~isempty(a), regexp(statMatrixColIDs,{[sprintf('_%s',bands{b}) '_Hilb']}));
            hilbVals{b} = statMatrix(:,hilbColLog);
            rmsVals{b} = zscore(conv(sqrt(conv(statMatrix(:,voltColLog).^2, ones(windowNdxs,1)/windowNdxs, 'same')), w, 'same'));
        else
            bands(b) = [];
        end
    end
    tetHilbVals{fl} = hilbVals;
    tetRmsVals{fl} = rmsVals;
    phasePowerNorm = cell(length(bands), length(bands));
    figure;
    spVals = nan(length(bands), length(bands));
    cLims = nan(length(bands),length(bands),2);
    for b1 = 1:length(bands)-1
        for b2 = 1:length(bands)
            if b1 >= b2
                continue
            else                
                phaseBins = linspace(-pi,pi,13);
                powerBins = linspace(-2,10,13);
                curPhasePowerVals = nan(length(powerBins)-1, length(phaseBins)-1);
                for pwr = 1:length(powerBins)-1
                    for ph = 1:length(phaseBins)-1
                        curPhasePowerLog = (rmsVals{b1}>=powerBins(pwr) & rmsVals{b1}<powerBins(pwr+1)) &...
                            (hilbVals{b1}>=phaseBins(ph) & hilbVals{b1}<phaseBins(ph+1));
                        curPhasePowerVals(pwr,ph) = mean(rmsVals{b2}(curPhasePowerLog));
                    end
                end
                phasePowerNorm{b1,b2} = curPhasePowerVals;
                spVals(b1,b2) = subplot(length(bands)-1, length(bands)-1, sub2ind([length(bands)-1 length(bands)-1], b1,b2-1));
                imagesc(phaseBins, powerBins, phasePowerNorm{b1,b2});
                set(spVals(b1,b2), 'ydir', 'normal');
%                 imagesc(phasePowerNorm{b1,b2});
                title(bands{b2});
                cLims(b1,b2,1) = min(get(spVals(b1,b2), 'clim'));
                cLims(b1,b2,2) = max(get(spVals(b1,b2), 'clim'));
                colorbar
                drawnow
            end
        end
    end
%     for r = 1:length(bands)
%         for c = 1:length(bands)-1
%             if r >= c
%                 continue;
%             else
%                 set(spVals(r,c), 'clim', [min(cLims(:)) max(cLims(:))]);
%             end
%         end
%     end
    colormap jet
    figTitle = annotation('textbox', 'position', [0.025 0.935 0.7 0.05], 'String', tetFiles{fl},...
        'linestyle', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
    set(gcf, 'PaperOrientation', 'landscape');
    print(gcf,'-fillpage');
end

%%
close all
load(fileNames{cellfun(@(a)~isempty(a), strfind(fileNames, 'BehaviorMatrix'))});

pokeInTimesMtx = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [0 0], 'PokeIn');
pokeInTimes = ExtractTrialData_SM(pokeInTimesMtx, behavMatrix(:,1));
pokeOutTimesMtx = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [0 0], 'PokeOut');
pokeOutTimes = ExtractTrialData_SM(pokeOutTimesMtx, behavMatrix(:,1));


                
    