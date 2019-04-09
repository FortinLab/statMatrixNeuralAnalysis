
% Phase-Power Coupling
%
%

%% Load Data
%% Extract Info and bands and shit idk
clear all
close all
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

for b = 1:length(bands)
    if sum(strcmp(fieldnames(windowSize), bands{b}))>=1
    else
        bands(b) = [];
    end
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
        windowNdxs = windowSize.(bands{b})*sampleRate;
        voltColLog = cellfun(@(a)~isempty(a), regexp(statMatrixColIDs,{[sprintf('_%s',bands{b}) '$']}));
        hilbColLog = cellfun(@(a)~isempty(a), regexp(statMatrixColIDs,{[sprintf('_%s',bands{b}) '_Hilb']}));
        hilbVals{b} = statMatrix(:,hilbColLog);
        rmsVals{b} = zscore(conv(sqrt(conv(statMatrix(:,voltColLog).^2, ones(windowNdxs,1)/windowNdxs, 'same')), w, 'same'));
    end
    tetHilbVals{fl} = hilbVals;
    tetRmsVals{fl} = rmsVals;
end
%%
for fl = 1:length(tetFiles)
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
%     print(gcf,'-fillpage');
end

%%
close all
load(fileNames{cellfun(@(a)~isempty(a), strfind(fileNames, 'BehaviorMatrix'))});
sampPsec = 1/mode(diff(behavMatrix(:,1)));

pokeInTimesMtx = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [0 0], 'PokeIn');
pokeInTimes = ExtractTrialData_SM(pokeInTimesMtx, behavMatrix(:,1));
pokeOutTimesMtx = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [0 0], 'PokeOut');
pokeOutTimes = ExtractTrialData_SM(pokeOutTimesMtx, behavMatrix(:,1));

odorVals = [pokeInTimesMtx.Odor];
posVals = [pokeInTimesMtx.Position];
inSeqLog = [pokeInTimesMtx.ItemItemDistance]==1;
performanceLog = [pokeInTimesMtx.Performance]==1;

timeWindows = [-0.5 0; 0 0.5];   
ndxRange = round(sampPsec*timeWindows);

trialPeriodIDs = [{'Pre-Trial'},{'Early Trial'},{'Late Trial'},{'Post Trial'}];

entrainerCol = 1;
phaseBins = linspace(-pi,pi,13);
powerBins = linspace(-3,10,13);
for tet = 1:length(tetHilbVals)
    figure;
    tempTrialHilbVals = cell(2,max(posVals),4);
    tempTrialRmsVals = cell(2,max(posVals),4);
    for trl = 1:length(pokeInTimes)
        pokeInNdx = find(behavMatrix(:,1)==pokeInTimes{trl});
        pokeOutNdx = find(behavMatrix(:,1)==pokeOutTimes{trl});
        if inSeqLog(trl) && performanceLog(trl)
            % Pre-Trial Period
            tempTrialHilbVals{1,posVals(trl),1} = [tempTrialHilbVals{1,posVals(trl),1}; cellfun(@(a)a(pokeInNdx+ndxRange(1,1):pokeInNdx+ndxRange(1,2)), tetHilbVals{tet}, 'uniformoutput', 0)];
            tempTrialRmsVals{1,posVals(trl),1} = [tempTrialRmsVals{1,posVals(trl),1}; cellfun(@(a)a(pokeInNdx+ndxRange(1,1):pokeInNdx+ndxRange(1,2)), tetRmsVals{tet}, 'uniformoutput', 0)];
            % Early-Trial Period
            tempTrialHilbVals{1,posVals(trl),2} = [tempTrialHilbVals{1,posVals(trl),2}; cellfun(@(a)a(pokeInNdx+ndxRange(2,1):pokeInNdx+ndxRange(2,2)), tetHilbVals{tet}, 'uniformoutput', 0)];
            tempTrialRmsVals{1,posVals(trl),2} = [tempTrialRmsVals{1,posVals(trl),2}; cellfun(@(a)a(pokeInNdx+ndxRange(2,1):pokeInNdx+ndxRange(2,2)), tetRmsVals{tet}, 'uniformoutput', 0)];
            % Late-Trial Period
            tempTrialHilbVals{1,posVals(trl),3} = [tempTrialHilbVals{1,posVals(trl),3}; cellfun(@(a)a(pokeOutNdx+ndxRange(1,1):pokeOutNdx+ndxRange(1,2)), tetHilbVals{tet}, 'uniformoutput', 0)];
            tempTrialRmsVals{1,posVals(trl),3} = [tempTrialRmsVals{1,posVals(trl),3}; cellfun(@(a)a(pokeOutNdx+ndxRange(1,1):pokeOutNdx+ndxRange(1,2)), tetRmsVals{tet}, 'uniformoutput', 0)];
            % Post-Trial Period
            tempTrialHilbVals{1,posVals(trl),4} = [tempTrialHilbVals{1,posVals(trl),4}; cellfun(@(a)a(pokeOutNdx+ndxRange(2,1):pokeOutNdx+ndxRange(2,2)), tetHilbVals{tet}, 'uniformoutput', 0)];
            tempTrialRmsVals{1,posVals(trl),4} = [tempTrialRmsVals{1,posVals(trl),4}; cellfun(@(a)a(pokeOutNdx+ndxRange(2,1):pokeOutNdx+ndxRange(2,2)), tetRmsVals{tet}, 'uniformoutput', 0)];
        else
        end
    end
    
    % Calculate CFC for all trials
    preTrialALLtrialsHilb = tempTrialHilbVals(1,:,1)'; preTrialALLtrialsHilb = preTrialALLtrialsHilb{:};
    preTrialALLtrials(:,:,1) = preTrialALLtrialsHilb;
    preTrialALLtrialsRMS = tempTrialRmsVals(1,:,1)'; preTrialALLtrialsRMS = preTrialALLtrialsRMS{:};
    preTrialALLtrials(:,:,2) = preTrialALLtrialsRMS;
    
    rlyTrialALLtrialsHilb = tempTrialHilbVals(1,:,2)'; rlyTrialALLtrialsHilb  = rlyTrialALLtrialsHilb{:};  
    rlyTrialALLtrials(:,:,1) = rlyTrialALLtrialsHilb;
    rlyTrialALLtrialsRMS = tempTrialRmsVals(1,:,2)'; rlyTrialALLtrialsRMS  = rlyTrialALLtrialsRMS{:};  
    rlyTrialALLtrials(:,:,2) = rlyTrialALLtrialsRMS;
        
    lteTrialALLtrialsHilb = tempTrialHilbVals(1,:,3)'; lteTrialALLtrialsHilb = lteTrialALLtrialsHilb{:};   
    lteTrialALLtrials(:,:,1) = lteTrialALLtrialsHilb;
    lteTrialALLtrialsRMS = tempTrialRmsVals(1,:,3)'; lteTrialALLtrialsRMS = lteTrialALLtrialsRMS{:}; 
    lteTrialALLtrials(:,:,2) = lteTrialALLtrialsRMS;
    
    postTrialALLtrialsHilb = tempTrialHilbVals(1,:,4)'; postTrialALLtrialsHilb = postTrialALLtrialsHilb{:};   
    postTrialALLtrials(:,:,1) = postTrialALLtrialsHilb;
    postTrialALLtrialsRMS = tempTrialRmsVals(1,:,4)'; postTrialALLtrialsRMS  = postTrialALLtrialsRMS{:}; 
    postTrialALLtrials(:,:,2) = postTrialALLtrialsRMS;
    
    trialPeriodData = [{cell2mat(preTrialALLtrials)}, {cell2mat(rlyTrialALLtrials)}, {cell2mat(lteTrialALLtrials)}, {cell2mat(postTrialALLtrials)}];
        
    for bnd = 1:length(bands)
        if bnd == entrainerCol
        else
            curPhasePowerVals = repmat({nan(length(powerBins)-1, length(phaseBins)-1)}, size(trialPeriodData));
            for pwr = 1:length(powerBins)-1
                for ph = 1:length(phaseBins)-1
                    for tprd = 1:length(trialPeriodData)
                        curPhasePowerLog = (trialPeriodData{tprd}(:,entrainerCol,2)>=powerBins(pwr) & trialPeriodData{tprd}(:,entrainerCol,2)<powerBins(pwr+1)) &...
                            (trialPeriodData{tprd}(:,entrainerCol,1)>=phaseBins(ph) & trialPeriodData{tprd}(:,entrainerCol,1)<phaseBins(ph+1));
                        curPhasePowerVals{tprd}(pwr,ph) = nanmean(trialPeriodData{tprd}(curPhasePowerLog,bnd,2));
%                         curPhasePowerVals{tprd}(pwr,ph) = nanmedian(trialPeriodData{tprd}(curPhasePowerLog,bnd,2));
                    end
                end
            end
            for tprd = 1:length(curPhasePowerVals)
                if entrainerCol==1 || bnd>entrainerCol
                    spVals(bnd,tprd) = subplot(length(bands)-1, length(curPhasePowerVals), sub2ind([length(bands)-1 length(curPhasePowerVals)], tprd,bnd-1));
                else
                    spVals(bnd,tprd) = subplot(length(bands)-1, length(curPhasePowerVals), sub2ind([length(bands)-1 length(curPhasePowerVals)], tprd,bnd));
                end
                imagesc(phaseBins, powerBins, curPhasePowerVals{tprd});
                set(spVals(bnd,tprd), 'ydir', 'normal');
                if tprd==1
                    ylabel(bands{bnd}, 'FontWeight', 'bold');
                end
                if bnd-entrainerCol ==1
                    title(trialPeriodIDs{tprd});
                end
                cLims(bnd,tprd,1) = min(get(spVals(bnd,tprd), 'clim')); 
                cLims(bnd,tprd,2) = max(get(spVals(bnd,tprd), 'clim')); 
                colorbar
                drawnow
            end
        end
    end
    
    
    for r = 1:length(bands)
        if r == entrainerCol
        else
            for c = 1:length(curPhasePowerVals)
                set(spVals(r,c), 'clim', [min(cLims(r,:)) max(cLims(r,:))]);
            end
        end
    end
    colormap jet
    figTitle = annotation('textbox', 'position', [0.025 0.935 0.7 0.05], 'String', tetFiles{tet},...
        'linestyle', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
    set(gcf, 'PaperOrientation', 'landscape');
    print(gcf,sprintf('%s_XfreqCouple', tetFiles{tet}(1:end-4)), '-dpdf', '-fillpage');
    close all
end
    