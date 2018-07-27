%% StatMatrixCreator
% Code to create the statMatrix data structure from recorded data
%% Identify experimental Type
% Tetrode Ensemble
% Multi-site LFP
% Multi-site Ensemble/LFP
expType = questdlg('What type of experiment are you compiling data from?',...
    'Identify Experimental Design',...
    'Single-site Ensemble',...
    'Multi-site LFP',...
    'ss Ensemble & ms LFP',...
    'Single-site Ensemble');
switch expType
    case 'Single-site Ensemble'
        exp = 1;
    case 'Multi-site LFP'
        exp = 2;
    case 'ss Ensemble & ms LFP'
        exp = 3;
    case ''
        disp('statMatrix creation cancelled');
end

%% Identify data source
% Open Ephys
% Plexon
%   -> Irvine
%   -> Boston
dataSource = questdlg('What system was used to collect the data?',...
    'Identify Data Acquisition System',...
    'Plexon',...
    'Open Ephys',...
    'Plexon');
switch dataSource
    case 'Plexon'
        plexonType = questdlg('Where was the data recorded?',...
            'Identify Plexon System',...
            'Irvine',...
            'Boston',...
            'Irvine');
        switch plexonType
            case 'Irvine'
                rig = 1;
            case 'Boston'
                plxDataCheck = questdlg('Make sure you have created and curated the plxData file for the recording session of interest',...
                    'plxData File Check',...
                    'Got it!',...
                    'Don''t got it, gotta make it',...
                    'Got it!');
                switch plxDataCheck
                    case 'Got it!'
                    case 'Don''t got it, gotta make it'
                        disp('statMatrix creation cancelled because of no plxData file');
                        return
                end
                rig = 2;
        end
        [fileName, filePath] = uigetfile('.plx','Identify .PLX Session File');
        if fileName == 0
            disp('No file selected');
            return
        end
        cd(filePath);
        plxFile = [filePath fileName];
        outputFileName = inputdlg('Determine File Suffix', 'Filename', 1, {fileName(1:end-4)});
    case 'Open Ephys'
        rig = 3;
        dir = uigetdir('Identify file directory for recording session');
        cd(dir);
        outputFileName = inputdlg('Determine File Suffix');
    case ''
        disp('statMatrix creation cancelled');
end

%% Data Type
% .plx Session File
% MountainSort cut .plx channels
%   -> with LFP in .plx session file
%   -> with LFP in .continuous OE files
dataType = questdlg('What type of data are you compiling?',...
    'Identify Data Storage',...
    'Plexon Session File',...
    'MountainSort Cut Files',...
    'Plexon Session File');
switch dataType
    case 'Plexon Session File'
        data = 1;
    case 'MountainSort Cut Files'
        questdlg({'Verify the cut, MS processed, .plx files have been named to match their channel names in the original plexon recordings.' '' 'If they don''t match, change them now'},...
            'File Name Verification', 'Ok', 'Ok');
        if rig == 1 || rig == 2
            data = 2;
        elseif rig == 3
            data = 3;
        end
    case ''
        disp('statMatrix creation cancelled');
end

%% Now create the statMatrix files
%% Create Behavior Matrix
if rig == 1                                                                 % Irvine .plx files
    [plxData] = SummarizePLXevents_SD(plxFile);
    behaviorData = plxData.Raw;
    summary = plxData.Summary;
    [behavMatrix, behavMatrixColIDs] = CreateBehaviorMatrix(rig, behaviorData, summary, outputFileName);
elseif rig == 2                                                             % Boston .plx files
    [fileName, filePath] = uigetfile('.mat','Select the plxData file for the recording session');
    load([filePath fileName]);
    trialType = questdlg('Which trials do you want to compile?',...
        'Determine trials to compile',...
        'All trials',...
        'Curated trials ONLY',...
        'Curated trials ONLY');
    switch trialType
        case 'All trials'
            behaviorData = plxData.Raw;
        case 'Curated trials ONLY'
            behaviorData = plxData.CuratedPLX;
    end
    summary = plxData.Summary;
    [behavMatrix, behavMatrixColIDs] = CreateBehaviorMatrix(rig, behaviorData, summary, outputFileName);
elseif rig == 3                                                             % Open Ephys files
    error('Open Ephys not implemented yet');
end
%% Create statMatrix
CreateNeuralMatrix(exp, data, rig, behavMatrix(:,1), summary, outputFileName)
    
%% Behavior Matrix Creation Functions
function [behavMatrix, behavMatrixColIDs] = CreateBehaviorMatrix(rig, behaviorData, summary, outputFileName)
    if rig == 1
        file = summary.PLXfile;
    elseif rig == 2
        file = summary.PlxFile;
    end

    chan = 1;
    [samp, ~, tetTS, fn, ~] = plx_ad_v(file, chan);
    while tetTS == -1
        chan = chan + 1;
        [samp, ~, tetTS, fn, ~] = plx_ad_v(file, chan);
    end

    % Create timestamp vector
    tsVect = ((0:(fn(1)-1))/samp)+tetTS(1);
    % Sometimes the file is composed of multiple fragments, in which case these
    % should be combined into a single file
    for fragNum = 2:length(tetTS)
        tsVect = [tsVect, ((0:(fn(fragNum)-1))/samp)+tetTS(fragNum)]; %#ok<AGROW>
    end
    tsVect(end+1) = tsVect(end)+(1/samp);

    seqLength = summary.SequenceLength;
    maxSeqLength = max([behaviorData.OrdinalPosition]);
    behVals = nan(length(tsVect)-1, seqLength + maxSeqLength + 5);
    behDataHeaders = cell(1,seqLength + maxSeqLength + 7);
    % Step through each sequence item/position and identify when odors were
    % presented
    for seq = 1:seqLength
        for pos = 1:maxSeqLength
            itemPresTimes = [behaviorData([behaviorData.SequenceItem]==seq).ItemPresentationTime];
            posPresTimes = [behaviorData([behaviorData.OrdinalPosition]==pos).ItemPresentationTime];
            behVals(:,seq) = histcounts(itemPresTimes, tsVect)';
            behVals(:,pos+seqLength) = histcounts(posPresTimes, tsVect)';
            behDataHeaders{seq} = ['Odor' num2str(seq)];
            behDataHeaders{pos+seqLength} = ['Position' num2str(pos)];
        end
    end
    inSeqOdorPres = [behaviorData([behaviorData.TranspositionDistance]==0).ItemPresentationTime];
    outSeqOdorPres = [behaviorData(~([behaviorData.TranspositionDistance]==0)).ItemPresentationTime];
    behVals(:,seqLength+maxSeqLength+1) = histcounts(inSeqOdorPres, tsVect)' - histcounts(outSeqOdorPres,tsVect)';
    behDataHeaders{maxSeqLength+seqLength+1} = 'InSeqLog';
    itmPresTimes = [behaviorData.ItemPresentationTime];
    trialPerformance = [behaviorData.Performance];
    corrTrials = itmPresTimes(logical(trialPerformance));
    corTrlHistCounts = histcounts(corrTrials, tsVect)';
    inCorrTrials = itmPresTimes(~logical(trialPerformance));
    inCorTrlHistCounts = histcounts(inCorrTrials, tsVect)';
    behVals(:,seqLength + maxSeqLength+2) = corTrlHistCounts + (inCorTrlHistCounts*-1);
    behDataHeaders{seqLength + maxSeqLength+2} = 'PerformanceLog';

    behVals(:,seqLength + maxSeqLength+3) = histcounts([behaviorData.OdorTrigPokeTime], tsVect)' - histcounts([behaviorData.OdorPokeWithdrawTime], tsVect)';
    behDataHeaders{seqLength + maxSeqLength+3} = 'PokeEvents';

    behVals(:,seqLength + maxSeqLength+4) = histcounts([behaviorData.FrontRewardTime], tsVect)';
    behDataHeaders{seqLength + maxSeqLength+4} = 'FrontReward';

    behVals(:,seqLength + maxSeqLength+5) = histcounts([behaviorData.BackRewardTime], tsVect)';
    behDataHeaders{seqLength + maxSeqLength+5} = 'BackReward';

    if isfield(behaviorData, 'ErrorSignalTime')
        behVals(:,seqLength + maxSeqLength+6) = histcounts([behaviorData.ErrorSignalTime], tsVect)';
        behDataHeaders{seqLength + maxSeqLength+6} = 'ErrorSignal';
    end

    [numChans, chanNames] = plx_event_names(file);
    findingStrobed = 1;
    while findingStrobed
        for chan = 1:numChans
            curChan = chanNames(chan,:);
            intVals = double(curChan);
            valLim = find(~(intVals==0), 1, 'last');
            strobedChanLog = strcmp(curChan(1:valLim), 'Strobed');
            if strobedChanLog
                [~, strobedTS, strobedSV] = plx_event_ts(file, curChan);
                [~, ~, ~, aniPosition] = plx_vt_interpret(strobedTS, strobedSV);
                aniPosition(aniPosition(:,1)<tsVect(1),:) = [];
                aniX = aniPosition(:,2);
                aniY = aniPosition(:,3);
                if size(aniPosition,2)>=5
                    for t = 1:size(aniPosition,1)
                        if (aniPosition(t,2)==0 && aniPosition(t,3)==0) &&...
                                (aniPosition(t,4)>0 && aniPosition(t,5)>0)
                            aniX(t) = aniPosition(t,4);
                            aniY(t) = aniPosition(t,5);
                        elseif (aniPosition(t,4)==0 && aniPosition(t,5)==0) &&...
                                (aniPosition(t,2)>0 && aniPosition(t,3)>0)
                            aniX(t) = aniPosition(t,2);
                            aniY(t) = aniPosition(t,3);
                        elseif (aniPosition(t,2)>0 && aniPosition(t,3)>0) &&...
                                (aniPosition(t,3)>0 && aniPosition(t,5)>0)
                            aniX(t) = mean([aniPosition(t,2) aniPosition(t,4)]);
                            aniY(t) = mean([aniPosition(t,3) aniPosition(t,5)]);
                        elseif (aniPosition(t,2)==0 && aniPosition(t,3)==0) &&...
                                (aniPosition(t,3)==0 && aniPosition(t,5)==0)
                            aniX(t) = 0;
                            aniY(t) = 0;
                        end
                    end
                end
                aniPosHistBins = find(histcounts(aniPosition(:,1), tsVect));
                if isfield(behaviorData, 'ErrorSignalTime')
                    behVals(aniPosHistBins,seqLength + maxSeqLength+7) = aniX';
                    behVals(aniPosHistBins,seqLength + maxSeqLength+8) = aniY';
                else
                    behVals(aniPosHistBins,seqLength + maxSeqLength+6) = aniX;
                    behVals(aniPosHistBins,seqLength + maxSeqLength+7) = aniY;
                end
                findingStrobed = 0;
            end
        end
    end
    if isfield(behaviorData, 'ErrorSignalTime')
        behDataHeaders{seqLength + maxSeqLength+7} = 'XvalRatMazePosition';
        behDataHeaders{seqLength + maxSeqLength+8} = 'YvalRatMazePosition';
    else
        behDataHeaders{seqLength + maxSeqLength+6} = 'XvalRatMazePosition';
        behDataHeaders{seqLength + maxSeqLength+7} = 'YvalRatMazePosition';
    end
    behavMatrix = [tsVect(1:end-1)', behVals];
    behavMatrixColIDs = [{'TimeBin'}, behDataHeaders];
    save([outputFileName{1} '_BehaviorMatrix.mat'], 'behavMatrix', 'behavMatrixColIDs');
    disp('Behavior data saved.');
end

%% Neural Matrix Creation Function
function CreateNeuralMatrix(exp, data, rig, tsVect, summary, outputFileName)
tsVect(end+1) = tsVect(end)+(mode(diff(tsVect)));
% Identify the source for the lfp Files 
if rig == 1 || rig == 2
    if rig == 1                                                         % Plexon collected in Irvine
        plxFile = summary.PLXfile;
    elseif rig == 2                                                     % Plexon collected in Boston
        plxFile = summary.PlxFile;
    end
    [~, plxADchanNames] = plx_adchan_names(plxFile);
    tetLFPchanNames = cell(size(plxADchanNames,1),1);
    for tet = 1:size(plxADchanNames,1)
        curADchan = plxADchanNames(tet,:);
        tetLFPchanNames{tet} = deblank(curADchan);
    end 
    [~, ~, ~, contCountFl] = plx_info(plxFile, 0);
    lfpDataLog = contCountFl ~= 0;
    tetLFPchanNames(~lfpDataLog) = [];
    tetLFPchanNames(cellfun(@(a)isempty(a), regexp(tetLFPchanNames, '^T([0-9]*)'))) = [];
elseif rig == 3                                                             % Open Ephys collection
    error('Open Ephys not implemented yet');
end
% Identify the Spiking Data
if exp == 1 || 3                                                            % Single-site ensemble data set
    if data == 1                                                            % Plexon session file
        % Unit and LFP data pulled from the same file
        [tsCountFl, ~, ~, ~] = plx_info(plxFile, 0);
        numUnis = sum(tsCountFl(:,2:end)~=0)-1;
        numUnis(numUnis==-1) = 0;
        [~,chanNames] = plx_chan_names(plxFile);
        tetNames = cell(size(chanNames,1),1);
        for t = 1:size(chanNames,1)
            curChanName = deblank(chanNames(t,:));
            tetNames{t} = curChanName(1:end-2);
        end
        tetsWithUnits = unique(tetNames(numUnis>=1));     
    elseif data == 2 || data == 3                                           % MountainSort cut files
        % Unit and LFP data pulled from different files
        %   Unit = individually cut .plx files
        %   LFP = PLX file (data==2) OR .CONTINUOUS files (data==3)
        spikeFileDir = uigetdir(cd, 'Where are the cut MountainSort .PLX files?');
        spikeFileDirFiles = dir([spikeFileDir '\']);
        spkFiles = {spikeFileDirFiles.name};
        spkFiles = spkFiles(logical(cellfun(@(a)~isempty(a), strfind(spkFiles, '.plx'))));
        if data == 2
            [tetNameStart, tetNameEnd] = regexp(spkFiles, '^T([0-9]*)');
        elseif data == 3
            % Determine naming convention for Open Ephys data files.
            error('Open Ephys not implemented yet');
        end            
        tetsWithUnits = cellfun(@(a,b,c)a(b:c), spkFiles, tetNameStart, tetNameEnd, 'uniformoutput', 0);
    end
end
ensembleMatrix = cell(size(tetsWithUnits));
ensembleMatrixColIDs = cell(size(tetsWithUnits));
ensembleUnitSummaries = cell(size(tetsWithUnits));
% Now run through each channel individually to identify things
% tetLFPchanNames
% tetsWithUnits
for chan = 1:length(tetLFPchanNames)
    if rig == 1 || rig == 2
        curADchan = tetLFPchanNames{chan};
        curTet = curADchan(1:end-2);
        
        fprintf('Found LFP data for %s on %s\n', curTet, curADchan);
        [samp, ~, ~, ~, tetV] = plx_ad_v(plxFile, curADchan);
        
        [tetRawHilb, tetRawLFP] = PhaseFreqDetectAbbr(tetV,tsVect(1:end-1));
        [tetThetaHilb, tetThetaLFP] = PhaseFreqDetectAbbr(tetV,tsVect(1:end-1),4,12);
        [tetLowBetaHilb, tetLowBetaLFP] = PhaseFreqDetectAbbr(tetV,tsVect(1:end-1),13,19);
        [tetBetaHilb, tetBetaLFP] = PhaseFreqDetectAbbr(tetV,tsVect(1:end-1),20,40);
        [tetLowGammaHilb, tetLowGammaLFP] = PhaseFreqDetectAbbr(tetV,tsVect(1:end-1),41,59);
        [tetHighGammaHilb, tetHighGammaLFP] = PhaseFreqDetectAbbr(tetV,tsVect(1:end-1),61,80);
        [tetRipHilb, tetRipLFP] = PhaseFreqDetectAbbr(tetV,tsVect(1:end-1),150,250);
        
        
        if ~isempty(find(strcmp(curTet, tetsWithUnits), 1))
            if data == 2 || data == 3
                tetFile = [spikeFileDir '\' spkFiles{find(strcmp(curTet, tetsWithUnits),1)}];
                [tsCountFl, ~, ~, ~] = plx_info(tetFile, 0);
                curNumUnis = sum(tsCountFl(:,2)>0)-1;            
                curChanNums = 1:4;
            else
                tetFile = plxFile;
                curUniLog = strcmp(tetNames, curTet);                
                curNumUnis = unique(numUnis(curUniLog));
                if length(curNumUnis)>1
                    error('Unit ID selection done plum goofed');
                end
                curChanNums = find(curUniLog);
            end
        else
            curNumUnis = 0;
        end
        
        % Now pull out spike data
        tetUnitTSs = cell(1,curNumUnis);
        tetUnitIDs = cell(1,curNumUnis);
        meanTemplates = repmat({cell(1,4)}, [1,curNumUnis]);
        meanTemplateValues = repmat({nan(3,2)}, [1, curNumUnis]);
        stdTemplates = repmat({cell(1,4)}, [1,curNumUnis]);
        meanSpkRt = cell(1,curNumUnis);
        spkWdth = cell(1,curNumUnis);
        phaseVals = struct('Mean', nan, 'Median', nan, 'R_Length', nan,...
            'AngDev', nan, 'CircStDev', nan, 'R_Test', [nan, nan]);
        phaseInfo = repmat({struct('Theta', phaseVals, 'LowBeta', phaseVals,...
            'Beta', phaseVals, 'LowGamma', phaseVals, 'HighGamma', phaseVals,...
            'Ripple', phaseVals)}, [1,curNumUnis]);
        for uni = 1:curNumUnis
            tetUnitIDs{uni} = [curTet '-U' num2str(uni)];
            fprintf('     Compiling unit %s\n', tetUnitIDs{uni});
            [~, tempTSs] = plx_ts(tetFile, curChanNums(1), uni);
            tetUnitTSs{uni} = histcounts(tempTSs, tsVect)';
            wire = 1;
            for tetChan = 1:length(curChanNums)
                [~, ~, ~, wave] = plx_waves_v(tetFile, curChanNums(tetChan), uni);
                meanTemplates{uni}{wire} = mean(wave);
                stdTemplates{uni}{wire} = std(wave);
                wire = wire+1;
            end
            % Calculate the spike width using the spike template
            % First find the largest template of the four wires
            templateHeights = cellfun(@(a)diff([min(a) max(a)]), meanTemplates{uni});
            largestTemplateNum = find(max(templateHeights)==templateHeights,1,'first');
            largestTemplate = meanTemplates{uni}{largestTemplateNum};
            meanTemplateValues{uni}(1,1) = largestTemplateNum;
            % Find the valley
            valleyNdx = find(largestTemplate==min(largestTemplate),1,'first');
            meanTemplateValues{uni}(2,:) = [valleyNdx, largestTemplate(valleyNdx)];
            % Bag the peak... following the valley
            peakNdx = find(largestTemplate(valleyNdx:end)==max(largestTemplate(valleyNdx:end)),1,'first') + (valleyNdx-1);
            meanTemplateValues{uni}(3,:) = [peakNdx, largestTemplate(peakNdx)];
            % Determine the width
            spkWdth{uni} = (peakNdx-valleyNdx)/samp;
            
            % Calculate the overall spike rate (spk/s)
            meanSpkRt{uni} = length(tempTSs)/(size(tetUnitTSs{uni},1)/samp);
            
            % Now pull out summary phase values
            % Theta
            thetaUniPhase = tetThetaHilb(tetUnitTSs{uni}==1);
            phaseInfo{uni}.Theta.Mean = circ_mean(thetaUniPhase);
            %         phaseInfo{uni}.Theta.Median = circ_median(thetaUniPhase);
            phaseInfo{uni}.Theta.R_Length = circ_r(thetaUniPhase);
            [phaseInfo{uni}.Theta.AngDev, phaseInfo{uni}.Theta.CircStDev] = circ_std(thetaUniPhase);
            [phaseInfo{uni}.Theta.R_Test(1), phaseInfo{uni}.Theta.R_Test(2)] = circ_rtest(thetaUniPhase);
            % Low Beta
            lowBetaUniPhase = tetLowBetaHilb(tetUnitTSs{uni}==1);
            phaseInfo{uni}.LowBeta.Mean = circ_mean(lowBetaUniPhase);
            %         phaseInfo{uni}.LowBeta.Median = circ_median(lowBetaUniPhase);
            phaseInfo{uni}.LowBeta.R_Length = circ_r(lowBetaUniPhase);
            [phaseInfo{uni}.LowBeta.AngDev, phaseInfo{uni}.LowBeta.CircStDev] = circ_std(lowBetaUniPhase);
            [phaseInfo{uni}.LowBeta.R_Test(1), phaseInfo{uni}.LowBeta.R_Test(2)] = circ_rtest(lowBetaUniPhase);
            % Beta
            betaUniPhase = tetBetaHilb(tetUnitTSs{uni}==1);
            phaseInfo{uni}.Beta.Mean = circ_mean(betaUniPhase);
            %         phaseInfo{uni}.Beta.Median = circ_median(betaUniPhase);
            phaseInfo{uni}.Beta.R_Length = circ_r(betaUniPhase);
            [phaseInfo{uni}.Beta.AngDev, phaseInfo{uni}.Beta.CircStDev] = circ_std(betaUniPhase);
            [phaseInfo{uni}.Beta.R_Test(1), phaseInfo{uni}.Beta.R_Test(2)] = circ_rtest(betaUniPhase);
            % Low Gamma
            lowGammaUniPhase = tetLowGammaHilb(tetUnitTSs{uni}==1);
            phaseInfo{uni}.LowGamma.Mean = circ_mean(lowGammaUniPhase);
            %         phaseInfo{uni}.LowGamma.Median = circ_median(lowGammaUniPhase);
            phaseInfo{uni}.LowGamma.R_Length = circ_r(lowGammaUniPhase);
            [phaseInfo{uni}.LowGamma.AngDev, phaseInfo{uni}.LowGamma.CircStDev] = circ_std(lowGammaUniPhase);
            [phaseInfo{uni}.LowGamma.R_Test(1), phaseInfo{uni}.LowGamma.R_Test(2)] = circ_rtest(lowGammaUniPhase);
            % High Gamma
            highGammaUniPhase = tetHighGammaHilb(tetUnitTSs{uni}==1);
            phaseInfo{uni}.HighGamma.Mean = circ_mean(highGammaUniPhase);
            %         phaseInfo{uni}.HighGamma.Median = circ_median(highGammaUniPhase);
            phaseInfo{uni}.HighGamma.R_Length = circ_r(highGammaUniPhase);
            [phaseInfo{uni}.HighGamma.AngDev, phaseInfo{uni}.HighGamma.CircStDev] = circ_std(highGammaUniPhase);
            [phaseInfo{uni}.HighGamma.R_Test(1), phaseInfo{uni}.HighGamma.R_Test(2)] = circ_rtest(highGammaUniPhase);
            % Ripple
            rippleUniPhase = tetRipHilb(tetUnitTSs{uni}==1);
            phaseInfo{uni}.Ripple.Mean = circ_mean(rippleUniPhase);
            %         phaseInfo{uni}.Ripple.Median = circ_median(rippleUniPhase);
            phaseInfo{uni}.Ripple.R_Length = circ_r(rippleUniPhase);
            [phaseInfo{uni}.Ripple.AngDev, phaseInfo{uni}.Ripple.CircStDev] = circ_std(rippleUniPhase);
            [phaseInfo{uni}.Ripple.R_Test(1), phaseInfo{uni}.Ripple.R_Test(2)] = circ_rtest(rippleUniPhase);
        end        
        unitSummary = struct('UnitName', tetUnitIDs, 'Spike_Features', meanTemplateValues,...
            'TemplateMean', meanTemplates, 'TemplateStDev', stdTemplates,...
            'Mean_SpikeRate', meanSpkRt, 'Spike_Width', spkWdth,...
            'Spike_Phase_Relations', phaseInfo); 
        statMatrix = [tsVect(1:end-1),tetRawLFP, tetRawHilb,...
            tetThetaLFP, tetThetaHilb,...
            tetLowBetaLFP, tetLowBetaHilb,...
            tetBetaLFP, tetBetaHilb,...
            tetLowGammaLFP, tetLowGammaHilb,...
            tetHighGammaLFP, tetHighGammaHilb,...
            tetRipLFP, tetRipHilb,...
            cell2mat(tetUnitTSs)]; %#ok<NASGU>
        statMatrixColIDs = [{'TimeBin'}, {[curTet '_LFP_Raw']}, {[curTet '_LFP_Raw_HilbVals']},...
            {[curTet '_LFP_Theta']}, {[curTet '_LFP_Theta_HilbVals']},...
            {[curTet '_LFP_LowBeta']}, {[curTet '_LFP_LowBeta_HilbVals']},...
            {[curTet '_LFP_Beta']}, {[curTet '_LFP_Beta_HilbVals']},...
            {[curTet '_LFP_LowGamma']}, {[curTet '_LFP_LowGamma_HilbVals']},...
            {[curTet '_LFP_HighGamma']}, {[curTet '_LFP_HighGamma_HilbVals']},...
            {[curTet '_LFP_Ripple']}, {[curTet '_LFP_Ripple_HilbVals']},...
            tetUnitIDs]; %#ok<NASGU>
        
        if curNumUnis >= 1
            ensembleMatrix{find(strcmp(curTet, tetsWithUnits), 1)} = cell2mat(tetUnitTSs);
            ensembleMatrixColIDs{find(strcmp(curTet, tetsWithUnits), 1)} = tetUnitIDs;
            ensembleUnitSummaries{find(strcmp(curTet, tetsWithUnits), 1)} = unitSummary;
        end
        
        save([outputFileName{1} '_' curTet '.mat'], 'statMatrix', 'statMatrixColIDs', 'unitSummary', '-v7.3');
        fprintf('          %s saved!\n', curTet);
    end
end

ensembleMatrix = [tsVect(1:end-1), cell2mat(ensembleMatrix)]; %#ok<NASGU>
ensembleMatrixColIDs = [{'TimeBin'}, ensembleMatrixColIDs{:}]; %#ok<NASGU>
ensembleUnitSummaries = cell2mat(ensembleUnitSummaries); %#ok<NASGU>
save([outputFileName{1} '_EnsembleMatrix.mat'], 'ensembleMatrix', 'ensembleMatrixColIDs', 'ensembleUnitSummaries', '-v7.3');
fprintf('%s ensembleMatrix saved!\n', outputFileName{1});
end