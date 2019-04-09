%% StatMatrixCreator
% Code to create the statMatrix data structure from recorded data
start = tic;
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
        recLoc = questdlg('What region is the ensemble recorded in?',...
            'Identify Recording Location',...
            'PFC',...
            'HC',...
            'Misc',...
            'PFC');
        switch recLoc
            case 'PFC'
                loc = 1;
            case 'HC'
                loc = 2;
            case 'Misc'
                loc = 3;
        end
    case 'Multi-site LFP'
        exp = 2;
    case 'ss Ensemble & ms LFP'
        exp = 3;
    case ''
        disp('statMatrix creation cancelled');
        return
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
                    case ''
                        disp('statMatrix creation cancelled');
                end
                rig = 2;
            case ''
                disp('statMatrix creation cancelled');
                return
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
        oeRig = questdlg('What setup was used to record the data?',...
            'Identify Recording Room',...
            '206',...
            '212',...
            '214',...
            '206');
        switch oeRig
            case '206'
                rig = 3;
            case '212'
                rig = 4;
            case '214'
                rig = 5;
            case ''
                disp('statMatrix creation cancelled');
                return
        end
        if exp == 2
            dir = uigetdir('Identify file directory for recording session');
            cd(dir);
            [fileName, filePath] = uigetfile('.mat', 'Identify the Channel ID Mapping file');
            chanMapFile = [filePath '\' fileName];
            load(chanMapFile);
        elseif exp == 3
            dir = uigetdir('Identify file directory for recording session');
            cd(dir);
        else
            error('Something''s wrong, start over (if you get this more than once talk to gabe');
        end
        outputFileName = inputdlg('Determine File Suffix', 'Filename', 1, {fileName(1:end-4)});
    case ''
        disp('statMatrix creation cancelled');
        return
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
    'Open Ephys .continuous',...
    'Plexon Session File');
switch dataType
    case 'Plexon Session File'
        data = 1;
    case 'MountainSort Cut Files'
        questdlg({'Verify the cut, MS processed, .plx files have been named to match their channel names in the original plexon recordings.' '' 'If they don''t match, change them now'},...
            'File Name Verification', 'Ok', 'Ok');
        if rig == 1 || rig == 2
            data = 2;
        elseif rig == 3 || rig == 4 || rig == 5
            data = 3;
        end
    case 'Open Ephys .continuous'
        data = 4;
    case ''
        disp('statMatrix creation cancelled');
        return
end
%% Create Log File
outfile = fopen([outputFileName{1} '_SMcreationLog.txt'], 'wt');
curTime = clock;
fprintf(outfile, 'Initialized %s at %i:%i.%i\n', date, curTime(4), curTime(5), round(curTime(6)));
fprintf(outfile, 'Working Directory = %s\n', cd);
if rig==1 || rig==2
    fprintf(outfile, '     Plexon File used = %s\n', plxFile);
elseif rig==3 || rig==4 || rig==5
    fprintf(outfile, '     OpenEphys Channel Map File Used = %s\n', chanMapFile);
end
fprintf(outfile, 'Experiment Type = %s\n', expType);
if exp==1
    fprintf(outfile, 'Recording Location = %s\n', recLoc);
end
fprintf(outfile, '     Data Source = %s\n', dataSource);
if rig==1 || rig==2
    fprintf(outfile, '     Rig = %s\n', plexonType);
elseif rig==3 || rig==4 || rig==5
    fprintf(outfile, '     Rig = %s\n', oeRig);
end
fprintf(outfile, '     Data Type = %s\n', dataType);
fprintf(outfile, '********************************************************\n');
%% Now create the statMatrix files
%% Create Behavior Matrix
if rig == 1                                                                 % Irvine .plx files
    [plxData] = SummarizePLXevents_SD(plxFile);
    behaviorData = plxData.Raw;
    summary = plxData.Summary;
    [behavMatrix, behavMatrixColIDs] = CreateBehaviorMatrixPLX(rig, behaviorData, summary, outputFileName, outfile);
elseif rig == 2                                                             % Boston .plx files
    [fileName, filePath] = uigetfile('.mat','Select the plxData file for the recording session');
    load([filePath fileName]);
    if fileName == 0
        fprintf(outfile, 'No file plxData file selected... statMatrix creator terminated.\n');
        error('No file selected, terminating statMatrix creation');
    else
        fprintf(outfile, 'Plexon summary file (plxData structure) used = %s\n', fileName);
        fprintf(outfile, '     File ID = %i\n', plxData.Summary.Identifier);
        fclose all;
        movefile([outputFileName{1} '_SMcreationLog.txt'], sprintf('%s_SMcreationLog_%i.txt', outputFileName{1}, plxData.Summary.Identifier));
        outfile = fopen(sprintf('%s_SMcreationLog_%i.txt', outputFileName{1}, plxData.Summary.Identifier), 'At');
    end        
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
    fprintf(outfile, '     Trial data used = %s\n', trialType);
    pathParts = strsplit(plxFile, '\');
    [behavMatrix, behavMatrixColIDs] = CreateBehaviorMatrixPLX(rig, behaviorData, summary, pathParts{end}, outputFileName, outfile);
elseif rig == 3 || rig == 4                                                         % Open Ephys files in a behavior rig
    [behavMatrix, behavMatrixColIDs] = CreateBehaviorMatrixOE(rig, outputFileName, outfile);
elseif rig == 5                                                                     % Open Ephys in main acoustic chamber (Tardis)
    error('No statMatrix creator code implemented for the tardis yet');
end
clc
%% Create statMatrix
if rig == 1 || rig == 2
    CreateNeuralMatrixPLX(exp, data, rig, loc, behavMatrix(:,1), summary, pathParts{end}, outputFileName, outfile)
else
    CreateNeuralMatrixOE(exp, data, rig, behavMatrix(:,1), chanMapStruct, outputFileName, outfile);
end

fprintf(outfile, 'StatMatrix creation complete. Process took %ih\n', round(toc(start),2)/3600);
fclose(outfile);
% SummarizeUnits_SM
%% Plexon Matrix Creation Functions
% Behavior Matrix
function [behavMatrix, behavMatrixColIDs] = CreateBehaviorMatrixPLX(rig, behaviorData, summary, fileName, outputFileName, outfile)
    global fileUseCheck % I know this is messy
    % Identify the PLX file being used for stuff
    if rig == 1
        file = summary.PLXfile;
    elseif rig == 2
        file = summary.PlxFile;
    end
    
    if strcmp(file, fileName)==0
        fileUseCheck = questdlg(sprintf('Filenames don''t match:\n plxFile = %s\n Selected = %s\n', file, fileName),...
            'Use selected files?',...
            'Yes',...
            'No',...
            'Yes');
        switch fileUseCheck
            case 'Yes'
                file = fileName;
                fprintf(outfile, 'plxData .plx file did not match selected file... used selected: %s\n', fileName);
            case 'No'
                fprintf(outfile, 'plxData .plx file did not match selected file... used plxData file: %s... code probably errored\n', file);
        end
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
    fprintf(outfile, 'Odor Counts\n');
    for seq = 1:seqLength
        itemPresTimes = [behaviorData([behaviorData.SequenceItem]==seq).ItemPresentationTime];
        behVals(:,seq) = histcounts(itemPresTimes, tsVect)';
        behDataHeaders{seq} = ['Odor' num2str(seq)];
        fprintf(outfile, '     Odor #%i = %i trials\n', seq, length(itemPresTimes));
    end
    fprintf(outfile, 'Position Counts\n');
    for pos = 1:maxSeqLength
        posPresTimes = [behaviorData([behaviorData.OrdinalPosition]==pos).ItemPresentationTime];
        behVals(:,pos+seqLength) = histcounts(posPresTimes, tsVect)';
        behDataHeaders{pos+seqLength} = ['Position' num2str(pos)];
        fprintf(outfile, '     Position #%i = %i trials\n', pos, length(posPresTimes));
    end
    inSeqOdorPres = [behaviorData([behaviorData.TranspositionDistance]==0).ItemPresentationTime];
    outSeqOdorPres = [behaviorData(~([behaviorData.TranspositionDistance]==0)).ItemPresentationTime];
    behVals(:,seqLength+maxSeqLength+1) = histcounts(inSeqOdorPres, tsVect)' - histcounts(outSeqOdorPres,tsVect)';
    behDataHeaders{maxSeqLength+seqLength+1} = 'InSeqLog';
    fprintf(outfile, '\nCompiling InSeq trials.....\n     %i trials were InSeq (%i%%)\n', length(inSeqOdorPres), round(length(inSeqOdorPres)/length(behaviorData),2)*100);
    itmPresTimes = [behaviorData.ItemPresentationTime];
    trialPerformance = [behaviorData.Performance];
    corrTrials = itmPresTimes(logical(trialPerformance));
    corTrlHistCounts = histcounts(corrTrials, tsVect)';
    inCorrTrials = itmPresTimes(~logical(trialPerformance));
    inCorTrlHistCounts = histcounts(inCorrTrials, tsVect)';
    behVals(:,seqLength + maxSeqLength+2) = corTrlHistCounts + (inCorTrlHistCounts*-1);
    behDataHeaders{seqLength + maxSeqLength+2} = 'PerformanceLog';
    fprintf(outfile, 'Compiling Performance.....\n     %i trials were correct (%i%%)\n', sum(trialPerformance), round(mean(trialPerformance),2)*100);

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
    fprintf(outfile, 'Behavior Matrix saved as %s_BehaviorMatrix.mat\n', outputFileName{1});
end

% Neural Matrix
function CreateNeuralMatrixPLX(exp, data, rig, loc, tsVect, summary, fileName, outputFileName, outfile)
    global fileUseCheck
    tsVect(end+1) = tsVect(end)+(mode(diff(tsVect)));
    % Identify the source for the lfp Files
    if rig == 1 || rig == 2
        if rig == 1                                                         % Plexon collected in Irvine
            plxFile = summary.PLXfile;
        elseif rig == 2                                                     % Plexon collected in Boston
            plxFile = summary.PlxFile;
        end
        switch fileUseCheck
            case 'Yes'
                plxFile = fileName;
        end
        [~, plxADchanNames] = plx_adchan_names(plxFile);
        tetLFPchanNames = cell(size(plxADchanNames,1),1);
        for tet = 1:size(plxADchanNames,1)
            curADchan = plxADchanNames(tet,:);
            tetLFPchanNames{tet} = deblank(curADchan);
        end
        [~, ~, ~, contCountFl] = plx_info(plxFile, 1);
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
            [tsCountFl, ~, ~, ~] = plx_info(plxFile, 1);
            numUnis = sum(tsCountFl(:,2:end)~=0)-1;
            numUnis(numUnis==-1) = 0;
            [~,chanNames] = plx_chan_names(plxFile);
            tetNames = cell(size(chanNames,1),1);
            for t = 1:size(chanNames,1)
                curChanName = deblank(chanNames(t,:));
                tetNames{t} = curChanName(1:end-2);
            end
            tetsWithUnits = unique(tetNames(numUnis>=1))';
        elseif data == 2 || data == 3                                           % MountainSort cut files
            % Unit and LFP data pulled from different files
            %   Unit = individually cut .plx files
            %   LFP = PLX file (data==2) OR .CONTINUOUS files (data==3)
            spikeFileDir = uigetdir(cd, 'Where are the cut MountainSort .PLX files?');
            fprintf(outfile, 'MountainSort pre-processed, manually curated .plx flies selected from: %s\n\n', [spikeFileDir '\']);
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
            fprintf(outfile, 'Processing Tetrode %s\n', curTet);
            fprintf(outfile, '     LFP data taken from %s\n', curADchan);

            [samp, ~, ~, ~, tetV] = plx_ad_v(plxFile, curADchan);
            
            if loc == 1
                bands = {'Raw', 'Theta', 'Alpha', 'Beta', 'LowGamma', 'HighGamma', 'Ripple'};
                bandLims = [[nan nan]; [4,8]; [9,12]; [16,32]; [35,55]; [61,80]; [150, 250]];
            elseif loc == 2 || loc == 3
                bands = {'Raw', 'Theta', 'LowBeta', 'Beta', 'LowGamma', 'HighGamma', 'Ripple'};
                bandLims = [[nan nan]; [4,12]; [13,19]; [20,40]; [40,59]; [61,80]; [150, 250]];
            end
            
            tetHilbVals = cell(1,length(bands));
            tetLFPvals = cell(1,length(bands));
            [tetHilbVals{1}, tetLFPvals{1}] = PhaseFreqDetectAbbr(tetV, tsVect(1:end-1));
            
            for b = 2:length(bands)
                [tetHilbVals{b}, tetLFPvals{b}] = PhaseFreqDetectAbbr(tetV, tsVect(1:end-1), bandLims(b,1), bandLims(b,2));
            end

            if ~isempty(find(strcmp(curTet, tetsWithUnits), 1))
                if data == 2 || data == 3
                    tetFile = [spikeFileDir '\' spkFiles{find(strcmp(curTet, tetsWithUnits),1)}];
                    [tsCountFl, ~, ~, ~] = plx_info(tetFile, 1);
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
                fprintf(outfile, '     Tetrode has %i units, extracted from %s\n', curNumUnis, tetFile);
            else
                curNumUnis = 0;
                fprintf(outfile, '     Tetrode has 0 units\n');
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
            phaseInfo = [];
            for b = 2:length(bands)
                phaseInfo.(bands{b}) = phaseVals;
            end
            phaseInfo = repmat({phaseInfo}, [1,curNumUnis]);
            for uni = 1:curNumUnis
                tetUnitIDs{uni} = [curTet '-U' num2str(uni)];
                fprintf('     Compiling unit %s\n', tetUnitIDs{uni});
                fprintf(outfile, '          Compiling %s.......', tetUnitIDs{uni});
                [~, tempTSs] = plx_ts(tetFile, curChanNums(1), uni);
                tetUnitTSs{uni} = histcounts(tempTSs, tsVect)';
                wire = 1;
                for tetChan = 1:length(curChanNums)
                    [~, ~, ~, wave] = plx_waves_v(tetFile, curChanNums(tetChan), uni);
                    meanTemplates{uni}{wire} = mean(wave);
                    stdTemplates{uni}{wire} = std(wave);
                    wire = wire+1;
                end
                fprintf(outfile, 'Waveform extraction complete, %i spikes found for %s\n', length(tempTSs), tetUnitIDs{uni});
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
                for b = 2:length(bands)
                    curUniPhase = tetHilbVals{b}(logical(tetUnitTSs{uni}));
                    phaseInfo{uni}.(bands{b}).Mean = circ_mean(curUniPhase);
%                     phaseInfo{uni}.(bands{b}).Median = circ_median(curUniPhase);
                    phaseInfo{uni}.(bands{b}).R_Length = circ_r(curUniPhase);
                    [phaseInfo{uni}.(bands{b}).AngDev, phaseInfo{uni}.(bands{b}).CircStDev] = circ_std(curUniPhase);
                    [phaseInfo{uni}.(bands{b}).R_Test(1), phaseInfo{uni}.(bands{b}).R_Test(2)] = circ_rtest(curUniPhase);
                end
            end
            unitSummary = struct('UnitName', tetUnitIDs, 'Spike_Features', meanTemplateValues,...
                'TemplateMean', meanTemplates, 'TemplateStDev', stdTemplates,...
                'Mean_SpikeRate', meanSpkRt, 'Spike_Width', spkWdth,...
                'LFPlimits', bandLims, 'Spike_Phase_Relations', phaseInfo);
            smBandMatrix = nan(length(tetLFPvals{1}), length(bands)*2);
            lfpColIDs = cell(1,length(bands)*2);
            bandNum = 0;
            for b = 1:2:length(bands)*2
                bandNum = bandNum+1;
                smBandMatrix(:,b) = tetLFPvals{bandNum};
                lfpColIDs{b} = sprintf('%s_LFP_%s', curTet, bands{bandNum});
                smBandMatrix(:,b+1) = tetHilbVals{bandNum};
                lfpColIDs{b+1} = sprintf('%s_LFP_%s_HilbVals', curTet, bands{bandNum});
            end
                
            statMatrix = [tsVect(1:end-1),...
                smBandMatrix...
                cell2mat(tetUnitTSs)]; %#ok<NASGU>
            statMatrixColIDs = [{'TimeBin'},...
                lfpColIDs,...
                tetUnitIDs]; %#ok<NASGU>

            if curNumUnis >= 1
                ensembleMatrix{find(strcmp(curTet, tetsWithUnits), 1)} = cell2mat(tetUnitTSs);
                ensembleMatrixColIDs{find(strcmp(curTet, tetsWithUnits), 1)} = tetUnitIDs;
                ensembleUnitSummaries{find(strcmp(curTet, tetsWithUnits), 1)} = unitSummary;
            end

            save([outputFileName{1} '_' curTet '_SM.mat'], 'statMatrix', 'statMatrixColIDs', 'unitSummary', '-v7.3');
            fprintf('          %s saved!\n', curTet);
            fprintf(outfile, '     %s saved as %s\n\n', curTet, sprintf('%s_%s.mat', outputFileName{1}, curTet));
        end
    end
    if ~isempty(tetsWithUnits(cellfun(@(a)isempty(a), ensembleMatrix)))
        lfpLessChan = tetsWithUnits{cellfun(@(a)isempty(a), ensembleMatrix)};
        fprintf(outfile, '%s was flagged as having units but has no LFP data. Channel skipped\n', lfpLessChan);
        ensembleUnitSummaries(cellfun(@(a)isempty(a), ensembleMatrix)) = [];
        ensembleMatrix(cellfun(@(a)isempty(a), ensembleMatrix)) = [];
    end
    if ~isempty(cell2mat(ensembleMatrix))
        ensembleMatrix = [tsVect(1:end-1), cell2mat(ensembleMatrix)]; %#ok<NASGU>
        ensembleMatrixColIDs = [{'TimeBin'}, ensembleMatrixColIDs{:}]; %#ok<NASGU>
        ensembleUnitSummaries = cell2mat(ensembleUnitSummaries); %#ok<NASGU>
        save([outputFileName{1} '_EnsembleMatrix.mat'], 'ensembleMatrix', 'ensembleMatrixColIDs', 'ensembleUnitSummaries', '-v7.3');
        fprintf('%s ensembleMatrix saved!\n', outputFileName{1});
        fprintf(outfile, '%s ensembleMatrix saved as %s\n', outputFileName{1}, sprintf('%s_EnsembleMatrix.mat', outputFileName{1}));
    end
end

%% Open Ephys Matrix Creation Functions
% Behavior Matrix
function [behavMatrix, behavMatrixColIDs] = CreateBehaviorMatrixOE(rig, outputFileName, outfile)
    % Event Channels:
    %       1) OdorA
    %       2) OdorB
    %       3) OdorC
    %       4) OdorD
    %       5) Front PB
    %       6) Front Water
    %       7) Rear Water
    %       8) Auditory Feedback
    % Pull filenames
    files = dir(cd);
    % Identify recording start time
    msgFlLog = cellfun(@(a)~isempty(a), strfind({files.name}, 'messages'));
    if sum(msgFlLog)==1
        msgFl = fopen(files(msgFlLog).name, 'r');
        txt = fgetl(msgFl);
        [recStartStart, recStartEnd] = regexp(txt, '^([0-9]*)');
        recStartTime = str2double(txt(recStartStart:recStartEnd));
        fclose(msgFl);
    elseif sum(msgFlLog)>=2
        error('Multiple ''messages'' files... why is it like this, it should not be like this!')
    else
        recStartTime = inputdlg('Specify recording start time');
    end
    
    % Find and load the ssnData file
    matFiles = {files(cellfun(@(a)~isempty(a), strfind({files.name}, '.mat'))).name};
    ssnDataFlLog = false(size(matFiles));
    for fl = 1:length(matFiles)
        variableInfo = who('-file', matFiles{fl});
        if sum(ismember(variableInfo, 'ssnData'))==1
            ssnDataFlLog(fl) = true;
        end
    end
    if sum(ssnDataFlLog)==1
        load(matFiles{ssnDataFlLog});
    elseif sum(ssnDataFlLog)>=2
        error('More than one .mat file.... this was prophesized! Tell gabe if you are not him!')
    else
        [ssnDataFileName,filePath] = uigetfile('.mat', 'Identify the ssnData file for this session');
        load([filePath ssnDataFileName]);
    end
    
    % Assess the ssnData structure; extract data and fill in missing info
    % (trial number and sequence number)
    if ~isfield(ssnData, 'TrialNum')
        trialNum = num2cell(1:length(ssnData));
    end
    if ~isfield(ssnData, 'SeqNum')
        seq = 0;
        seqNum = cell(size(ssnData));
        for trl = 1:length(ssnData)
            if ssnData(trl).TrialPosition == 1
                seq = seq+1;
            end
            seqNum{trl} = seq;
        end
    end
    sequenceLength = ssnData(1).Settings.SequenceLength;
    bufferPrd = ssnData(1).Settings.ShortPokeBuffer;
    if isfield(ssnData(1).Settings, 'ShortPokeBufferDur')
        bufferDur = ssnData(1).Settings.ShortPokeBufferDur;
    else
        bufferDur = 0.3;
    end
        
    adcFiles = {files(cellfun(@(a)~isempty(a),regexp({files.name}, '100_ADC([1-9]*)'))).name};
    adcEventTimes = cell(length(adcFiles),2);
    for fl = 1:length(adcFiles)
        fprintf('Pulling event record from %s...', adcFiles{fl});
        if fl==1
            [tempContData,~,info] = load_open_ephys_data_faster(adcFiles{1});
            sampleRate = info.header.sampleRate;
        else
            [tempContData,~,~] = load_open_ephys_data_faster(adcFiles{fl});
        end
%         if rig == 3
            [~, offTimes] = findpeaks([0; diff(tempContData)],sampleRate, 'MinPeakProminence', 0.25);
            adcEventTimes{fl,2} = offTimes + recStartTime;
            [~, onTimes] = findpeaks([0; diff(tempContData)*-1],sampleRate, 'MinPeakProminence', 0.25);
            adcEventTimes{fl,1} = onTimes + recStartTime;
%         elseif rig == 4
%             [~, adcEventTimes{fl,1}] = findpeaks([0; diff(tempContData)],sampleRate, 'MinPeakProminence', 0.25);
%             [~, adcEventTimes{fl,2}] = findpeaks([0; diff(tempContData)*-1],sampleRate, 'MinPeakProminence', 0.25);
%         end
        if length(adcEventTimes{fl,1}) + 1 == length(adcEventTimes{fl,2})
            adcEventTimes{fl,2}(end) = [];
        elseif length(adcEventTimes{fl,2}) + 1 == length(adcEventTimes{fl,1})
            adcEventTimes{fl,1}(end) = [];
        end
        fprintf('timestamps found\n');
    end   
    odorOnTimes = adcEventTimes(1:4,1);
    odorOnMatrix = cellfun(@(a,b) [a(:), ones(length(a),1)*b], odorOnTimes, num2cell(1:4)', 'uniformoutput', 0);
    odorOnMatrix = sortrows(cell2mat(odorOnMatrix));
    
    % Check the number of odor events vs number of trials in the ssnData
    % structure.
    if length(ssnData) ~= size(odorOnMatrix,1)
        if sum(odorOnMatrix(end-3:end,2) == [1:4]')==4 && sum(diff(odorOnMatrix(end-3:end,1))) == 0 % Captured odor pressure release
            odorOnMatrix = odorOnMatrix(1:end-4,:);
            if length(ssnData) ~= size(odorOnMatrix,1)
                error('Something no correct... odors no match!')
            elseif sum([ssnData.Odor] - odorOnMatrix(:,2)') ~=0
            end
        else
            error('Odor lengths do not match... something is wrong');
        end
    end                
    odorOnMatrix(:,3) = [ssnData.TrialPosition]';
    
    % Now create the timestamp vector used to bin the events.
    tsVect = recStartTime:1/sampleRate:recStartTime+(length(tempContData)/sampleRate);
    % If the sampleRate was higher than 1k (which is an unofficial standard
    % given the previously collected data). The timestamp vector should be
    % downsampled to create what a 1k sample rate would have seen.
    %       NOTE: The ADC channels should NOT be downsampled, as doing so
    %           may cause issues 
    if sampleRate ~= 1000
        dsRate = sampleRate/1000;
        tsVect = downsample(tsVect,dsRate);
    end
        
    % Create the behavMatrix
    behavMatrix = [tsVect(1:end-1)',...
        histcounts(odorOnMatrix(odorOnMatrix(:,2)==1,1), tsVect)',...
        histcounts(odorOnMatrix(odorOnMatrix(:,2)==2,1), tsVect)',...
        histcounts(odorOnMatrix(odorOnMatrix(:,2)==3,1), tsVect)',...
        histcounts(odorOnMatrix(odorOnMatrix(:,2)==4,1), tsVect)',...
        histcounts(odorOnMatrix(odorOnMatrix(:,3)==1,1), tsVect)',...
        histcounts(odorOnMatrix(odorOnMatrix(:,3)==2,1), tsVect)',...
        histcounts(odorOnMatrix(odorOnMatrix(:,3)==3,1), tsVect)',...
        histcounts(odorOnMatrix(odorOnMatrix(:,3)==4,1), tsVect)',...
        zeros(length(tsVect)-1,6)];
    behavMatrixColIDs = [{'TimeBin'}, {'Odor1'}, {'Odor2'}, {'Odor3'}, {'Odor4'},...
        {'Position1'}, {'Position2'}, {'Position3'}, {'Position4'},...
        {'InSeqLog'}, {'PerformanceLog'}, {'PokeEvents'}, {'FrontReward'},...
        {'BackReward'}, {'ErrorSignal'}];
    
    % This may seem odd to do things this way (below) but it allows
    % flexibility in the ColIDs above, without having to change all the
    % tons of stuff below to make the statMatrix with it.
    isCol = strcmp(behavMatrixColIDs, 'InSeqLog');
    perfCol = strcmp(behavMatrixColIDs, 'PerformanceLog');
    pokeCol = strcmp(behavMatrixColIDs, 'PokeEvents');
    fRwdCol = strcmp(behavMatrixColIDs, 'FrontReward');
    rRwdCol = strcmp(behavMatrixColIDs, 'BackReward');
    errCol = strcmp(behavMatrixColIDs, 'ErrorSignal');
    
    % Fill In InSeqLog Column
    isLog = [ssnData.TranspositionDistance]==0;
    behavMatrix(:,isCol) = histcounts(odorOnMatrix(isLog,1), tsVect);
    behavMatrix(:,isCol) =  behavMatrix(:,isCol) - histcounts(odorOnMatrix(~isLog,1), tsVect)';
    
    % Fill In PerformanceLog Column
    perfLog = [ssnData.Performance]==1;
    behavMatrix(:,perfCol) = histcounts(odorOnMatrix(perfLog,1), tsVect);
    behavMatrix(:,perfCol) = behavMatrix(:,perfCol) - histcounts(odorOnMatrix(~perfLog,1), tsVect)';
    
    % Fill In Front Reward Column
    behavMatrix(:,fRwdCol) = histcounts(adcEventTimes{6,1}, tsVect);
    behavMatrix(:,fRwdCol) = behavMatrix(:,fRwdCol) - histcounts(adcEventTimes{6,2}, tsVect)';
    
    % Fill in Rear Reward Column
    behavMatrix(:,rRwdCol) = histcounts(adcEventTimes{7,1}, tsVect);
    behavMatrix(:,rRwdCol) = behavMatrix(:,rRwdCol) - histcounts(adcEventTimes{7,2}, tsVect)';
    
    % Extract Auditory Signal Data (Buzzer... error signal & trial/sequence
    % available signals based on temporal characteristics) & Fill in Error
    % Signal column
    audSigDurs = diff(cell2mat(adcEventTimes(8,:)),1,2);
    behavMatrix(:,errCol) = histcounts(adcEventTimes{8,1}(audSigDurs>0.15), tsVect);
    behavMatrix(:,errCol) = behavMatrix(:,errCol) - histcounts(adcEventTimes{8,2}(audSigDurs>0.15), tsVect)';
    
    % Extract Poke Events
    pokeTSs = [adcEventTimes{5,:}];   
    % Calculate Inter-Poke-Intervals
    ipi = nan(size(pokeTSs,1),1);
    for poke = 1:length(ipi)-1
        ipi(poke) = pokeTSs(poke+1,1) - pokeTSs(poke,2);
    end
    pokeInTSs = nan(size(ssnData));
    pokeOutTSs = nan(size(ssnData));
    for trl = 1:size(odorOnMatrix,1)
        curOdorTS = odorOnMatrix(trl,1);
        pokeNdx = find(pokeTSs(:,1)<curOdorTS,1,'last');
        pokeInTSs(trl) = pokeTSs(pokeNdx,1);
        
        if trl == size(odorOnMatrix,1)
            trialPokes = pokeTSs(pokeTSs(:,1)>=pokeTSs(pokeNdx,1),:);
            trialIPI = ipi(pokeTSs(:,1)>=pokeTSs(pokeNdx,1),:);
        else
            trialPokes = pokeTSs(pokeTSs(:,1)>=pokeTSs(pokeNdx,1) & pokeTSs(:,1)<odorOnMatrix(trl+1,1),:);
            trialIPI = ipi(pokeTSs(:,1)>=pokeTSs(pokeNdx,1) & pokeTSs(:,1)<odorOnMatrix(trl+1,1),:);
            trialPokes(end,:) = []; % Remove the last poke from here because that's what initiates the next trial.
            trialIPI(end,:) = [];
        end
        
        trialPokesMatrix = [trialPokes, trialPokes-trialPokes(1), trialIPI];
                
        if (round(trialPokesMatrix(1,4),4) == round(ssnData(trl).PokeDuration,4) ||...  % IF the poke durations match
                round(trialPokesMatrix(1,4),5) == round(ssnData(trl).PokeDuration,5))          %... once floating point values are removed.
            pokeOutTSs(trl) = pokeTSs(pokeNdx,2);                                   % THEN use the poke out value associated with that poke in! ezpz
        else                                                                    % HOWEVER, if they don't match, check the trial based on it's type
            if ssnData(trl).TranspositionDistance == 0                              % IF the trial is InSeq
                if ssnData(trl).Performance == 1                                        % IF they got it correct
                    if (trialPokesMatrix(1,4) > ssnData(trl).TargetPokeDur)                % IF they held for longer than the trial's target duration
                        pokeOutTSs(trl) = pokeTSs(pokeNdx,2);                                       % THEN use the poke out value associated with that poke in!
                    elseif trialPokesMatrix(1,4) < ssnData(trl).TargetPokeDur              % ELSEIF they didn't hold long enough on the poke that started the trial... they should have triggerd the buffer
                        if size(trialPokesMatrix,1)>=2 &&...                                         % IF there was more than one poke on that trial (which there SHOULD be)
                                trialPokesMatrix(1,4) < bufferPrd &&...                                         % AND they withdrew during the buffer period
                                trialPokesMatrix(1,5) < bufferDur                                               % AND they poked back in within the buffer duration
                            curPokeDur = trialPokesMatrix(1,4);                                             % THEN step through the poke matrix to identify when they crossed threshold!
                            pk = 1;
                            while curPokeDur < ssnData(trl).TargetPokeDur &&...
                                    pk<=size(trialPokesMatrix,1)
                                pk = pk + 1;
                                curPokeDur = trialPokesMatrix(pk,4);
                            end
                            % After stepping through it SHOULD be at a point
                            % where the hold duration is correctly beyond the
                            % target.
                            if curPokeDur > ssnData(trl).TargetPokeDur
                                pokeOutTSs(trl) = trialPokesMatrix(pk,2);
                            else
                                % If it isn't... something's wrong.
                                error('Code no should go here... something wrong');
                            end
                        elseif (round(trialPokesMatrix(1,4),2) >= round(ssnData(trl).TargetPokeDur,2))
                            pokeOutTSs(trl) = pokeTSs(pokeNdx,2);
                        else
                            error('Code no should go here also... something wrong');
                        end
                    else
                        error('Code no should go here either... something wrong');
                    end
                elseif ssnData(trl).Performance == 0
                    % With InSeq trials, the error signal should happen
                    % AFTER they withdraw
                    audSigTS = adcEventTimes{8,1}(find(adcEventTimes{8,1}>trialPokesMatrix(1,1), 1, 'first'));
                    pokeOutTSs(trl) = trialPokesMatrix(find(trialPokesMatrix(:,2)<audSigTS,1,'last'),2);                    
                else
                    error('Huh?')
                end
            else
                if size(trialPokesMatrix,1) == 1                            % IF there's only one poke...
                    pokeOutTSs(trl) = trialPokesMatrix(1,2);                    % THEN just use the poke out value because that's all there is! ezpz
                else                                                        % ELSE there must be more than one poke... so time to examine the trial events
                    if (ssnData(trl).Performance == 1)                          % IF they got the trial correct
                        nxtRwdTS = adcEventTimes{6,1}(find(adcEventTimes{6,1}>trialPokesMatrix(1,1),1,'first'));
                        if trl~=size(odorOnMatrix,1) &&...                  % IF it's not the last trial
                                nxtRwdTS < odorOnMatrix(trl+1,1)                % AND the next reward timestamp comes before the next odor presentation (which it should in this case)
                            % Remove any poke events that happened after
                            % the reward ... I am confident in the logic of
                            % my control code that there shouldn't be any
                            trialPokesMatrix(trialPokesMatrix(:,1)>nxtRwdTS,:) = [];
                            % Then use the final poke out value from the
                            % edited trial pokes matrix as the poke out
                            % time
                            pokeOutTSs(trl) = trialPokesMatrix(end,2);
                        else
                            error('Code should not go here.');
                        end
                    elseif ssnData(trl).Performance == 0
                        % If they got it wrong, find the first triggering
                        % of the error channel after the poke initiation
                        audSigTS = adcEventTimes{8,1}(find(adcEventTimes{8,1}>trialPokesMatrix(1,1), 1, 'first'));
                        % With OutSeq trials the error should trigger after
                        % the decision threshold passes, i.e. poke
                        % withdrawal should follow the error signal trigger
                        pokeOutTSs(trl) = trialPokesMatrix(find(trialPokesMatrix(:,2)>audSigTS,1,'first'),2);
                    else
                        error('Impossible! Inconceivable! Unexpected that the code would go here!');
                    end
                end
            end
        end
    end
    % Now fill in the Poke events column
    behavMatrix(:,pokeCol) = histcounts(pokeInTSs, tsVect);
    behavMatrix(:,pokeCol) = behavMatrix(:,pokeCol) - histcounts(pokeOutTSs, tsVect)';  
    
    % Now Save it all!    
    save([outputFileName{1} '_BehaviorMatrix.mat'], 'behavMatrix', 'behavMatrixColIDs');
    disp('Behavior data saved.');
    fprintf(outfile, 'Behavior Matrix saved as %s_BehaviorMatrix.mat\n', outputFileName{1});
end
% 
function CreateNeuralMatrixOE(exp, data, rig, tsVect, chanMapStruct, outputFileName, outfile)
    % Define the LFP bands to be used in the statMatrix
    lfpBands = [{'Raw'}, {[]},...
        {'Theta'}, {[4 12]},...
        {'LowBeta'}, {[13 19]},...
        {'Beta'}, {[20 40]},...
        {'LowGamma'}, {[41 59]},...
        {'HighGamma'}, {[61 80]},...
        {'Ripple'}, {[150 250]}];
    
    tsVect(end+1) = tsVect(end)+(mode(diff(tsVect)));
    % For this one... right now focus on the multi-site LFP code... put
    % errors in for other exp codes.
    if exp == 2 
    elseif exp == 3 || data == 3
        error('Code not implemented yet to compile ensemble and multi-site LFP data');
    else
        error('Incorrect experiment selection... start over');
    end
    if rig == 3 || rig == 4 || rig == 5
    else
        error('Incorrect rig selection... start over');
    end
    % Determine regional groupings
    fileIDs = chanMapStruct.FileIDs;
    
    for fl = 1:size(fileIDs,1)
        [tempContData,~,info] = load_open_ephys_data_faster(fileIDs{fl,1});
        sampleRate = info.header.sampleRate;
        if sampleRate ~= 1000
            dsRate = sampleRate/1000;
            tempContData = downsample(tempContData,dsRate);
            if length(tempContData) == length(tsVect)
                tempContData = tempContData(1:end-1);
            elseif length(tempContData) ~= length(tsVect)-1
                error('Something wrong!');
            end
        end
        statMatrix = nan(length(tsVect)-1, length(lfpBands));
        statMatrixColIDs = cell(1,length(lfpBands));
        for bnd = 2:2:length(lfpBands)
            if isempty(lfpBands{bnd})
                [statMatrix(:,bnd), statMatrix(:,bnd-1)] = PhaseFreqDetectAbbr(tempContData,tsVect(1:end-1));
            else
                [statMatrix(:,bnd), statMatrix(:,bnd-1)] = PhaseFreqDetectAbbr(tempContData,tsVect(1:end-1),lfpBands{bnd}(1),lfpBands{bnd}(2));
            end
            statMatrixColIDs{bnd-1} = [fileIDs{fl,2} '_LFP_' lfpBands{bnd-1}];
            statMatrixColIDs{bnd} = [fileIDs{fl,2} '_LFP_' lfpBands{bnd-1} '_HilbVals'];
        end
        statMatrix = [tsVect(1:end-1),statMatrix];
        statMatrixColIDs = [{'TimeBin'},statMatrixColIDs];
        
        save([outputFileName{1} '_' fileIDs{fl,2} '_SM.mat'], 'statMatrix', 'statMatrixColIDs','-v7.3');
        fprintf('          %s saved!\n', fileIDs{fl,2});
        fprintf(outfile, '     %s saved as %s\n\n', fileIDs{fl,2}, sprintf('%s_%s_SM.mat', outputFileName{1}, fileIDs{fl,2}));
    end
end