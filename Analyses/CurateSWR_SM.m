aniIDs = [{'Barat'},...
    {'Buchanan'},...
    {'Mitt'},...
    {'Stella'},...
    {'SuperChris'}];
dataDir = 'D:\WorkBigDataFiles\CA1 Data\1. WellTrained session\';
aniDirs = cellfun(@(a)sprintf('%s%s\\SWR Tets\\',dataDir, a), aniIDs, 'uniformoutput', 0);
swChans = [{16},...
    {21},...
    {2},...
    {18},...
    {12}];
ripChans = [{21},...
    {18},...
    {18},...
    {14},...
    {15}];
aniInfo = struct('ID', aniIDs, 'Directory', aniDirs, 'SWchan', swChans, 'RIPchan', ripChans,...
    'SWepocs', [], 'SWtrace', [], 'SWpower', [],...
    'RIPepocs', [], 'RIPtrace', [], 'RIPpower', [],...
    'SWRepocs', []);

swThresh = [0 3];
ripThresh = [0 3];
mergeThresh = 15;

odorColors = [44/255, 168/255, 224/255;...
    154/255, 133/255, 122/255;...
    9/255, 161/255, 74/255;...
    128/255, 66/255, 151/255;...
    241/255, 103/255, 36/255];

%%
for a = 1:length(aniInfo)
    %% Load Files & Load Behavior/Ensemble Data
    aniFiles = dir(aniInfo(a).Directory);
    behMatFile = aniFiles(cellfun(@(a)~isempty(a), strfind({aniFiles.name}, 'BehaviorMatrix'))).name;
    behav = load([aniInfo(a).Directory behMatFile]);
    behavMat = OrganizeTrialData_SM(behav.behavMatrix, behav.behavMatrixColIDs, [0 0], 'PokeIn');
    sessionMatrix = [[behavMat.PokeInIndex]', [behavMat.PokeOutIndex]', [behavMat.Performance]', [behavMat.Position]', [behavMat.Odor]'];
    clear behav behavMat
    nsmblMatFile = aniFiles(cellfun(@(a)~isempty(a), strfind({aniFiles.name}, 'EnsembleMatrix'))).name;
    ensemble = load([aniInfo(a).Directory nsmblMatFile]);
    
    %% Evaluate Sharpwave channel & events
    swFile = aniFiles(cell2mat(cellfun(@(a)~isempty(a),regexp({aniFiles.name}, ['\w*T' mat2str(aniInfo(a).SWchan) '_\w*SM.mat\>']),'uniformoutput', 0))).name;
    sw = load([aniInfo(a).Directory swFile]);
    [swEpocs, swLFP, swPOW] = SharpwaveDetection(sw.statMatrix(:,2), 1/mode(diff(sw.statMatrix(:,1))), swThresh, mergeThresh);
    swLog = false(size(sw.statMatrix(:,1)));
    for sws = 1:size(swEpocs,1)
        swLog(swEpocs(sws,1):swEpocs(sws,2)) = true;
    end
    swTrace = nan(size(sw.statMatrix(:,1)));
    swTrace(swLog) = sw.statMatrix(swLog,2);
    aniInfo(a).SWepocs = swEpocs;
    aniInfo(a).SWtrace = swLFP;
    aniInfo(a).SWpower = swPOW;
    
    swPV = nan(size(ensemble.ensembleMatrix,2)-1, size(swEpocs,1));
    for sws = 1:size(swEpocs,1)
        swPV(:,sws) = mean(ensemble.ensembleMatrix(swEpocs(sws,1):swEpocs(sws,2),2:end),1);
    end
        
    %% Evaluate Ripple channel & events
    ripFile = aniFiles(cell2mat(cellfun(@(a)~isempty(a),regexp({aniFiles.name}, ['\w*T' mat2str(aniInfo(a).RIPchan) '_\w*SM.mat\>']),'uniformoutput', 0))).name;
    rip = load([aniInfo(a).Directory ripFile]);
    [ripEpocs, ripLFP, ripPOW] = RippleDetection(rip.statMatrix(:,2), 1/mode(diff(rip.statMatrix(:,1))), ripThresh, mergeThresh);
    ripLog = false(size(rip.statMatrix(:,1)));
    for rips = 1:size(ripEpocs,1)
        ripLog(ripEpocs(rips,1):ripEpocs(rips,2)) = true;
    end
    ripTrace = nan(size(rip.statMatrix(:,1)));
    ripTrace(ripLog) = rip.statMatrix(ripLog,2);
    aniInfo(a).RIPepocs = ripEpocs;
    aniInfo(a).RIPtrace = ripLFP;
    aniInfo(a).RIPpower = ripPOW;
    
    ripPV = nan(size(ensemble.ensembleMatrix,2)-1, size(ripEpocs,1));
    for rips = 1:size(ripEpocs,1)
        ripPV(:,rips) = mean(ensemble.ensembleMatrix(ripEpocs(rips,1):ripEpocs(rips,2),2:end),1);
    end
    
    %% Identify Sharpwave Ripple events
    swrLog = swLog + ripLog;
    swrStart = find(diff(swrLog==2)==1);
    swrEnd = find(diff(swrLog==2)==-1);
    swrWindows = [swrStart nan(length(swrStart),1)];
    for swr = 1:length(swrStart)
        tempSWRstart = swrStart(swr);
        swrWindows(swr,1) = min([swEpocs(tempSWRstart>=swEpocs(:,1) & tempSWRstart<=swEpocs(:,2),1), ripEpocs(tempSWRstart>=ripEpocs(:,1) & tempSWRstart<=ripEpocs(:,2),1), tempSWRstart]);
        
        tempSWRend = swrEnd(find(swrEnd>tempSWRstart,1,'first'));
        swrWindows(swr,2) = max([swEpocs(tempSWRend>=swEpocs(:,1) & tempSWRend<=swEpocs(:,2),2), ripEpocs(tempSWRend>=ripEpocs(:,1) & tempSWRend<=ripEpocs(:,2),2), tempSWRend]);
    end
    for swr = 2:size(swrWindows,1)
        if ~isnan(swrWindows(swr,1)) && ~isnan(swrWindows(swr-1,1))
            if swrWindows(swr,1) - swrWindows(swr-1,2) <= mergeThresh
                swrWindows(swr,1) = swrWindows(swr-1,1);
                swrWindows(swr-1,:) = nan;
            end
        end
    end
    swrWindows(isnan(swrWindows(:,1)),:) = [];
    aniInfo(a).SWRepocs = swrWindows;
    
    newSWtrace = nan(size(sw.statMatrix(:,1)));
    newRIPtrace = nan(size(rip.statMatrix(:,1)));
    for swr = 1:size(swrWindows,1)
        newSWtrace(swrWindows(swr,1):swrWindows(swr,2)) = sw.statMatrix(swrWindows(swr,1):swrWindows(swr,2),2);
        newRIPtrace(swrWindows(swr,1):swrWindows(swr,2)) = rip.statMatrix(swrWindows(swr,1):swrWindows(swr,2),2);
    end
    
    swrPV = nan(size(ensemble.ensembleMatrix,2)-1, size(swrWindows,1));
    swrActv = nan(1,size(swrWindows,1));
    for swrs = 1:size(swrWindows,1)
        swrPV(:,swrs) = mean(ensemble.ensembleMatrix(swrWindows(swrs,1):swrWindows(swrs,2),2:end),1);
        swrActv(swrs) = mean(swrPV(:,swrs)~=0);
    end
    
    %% Summarize SWR Events
    % Tabulate SWR vs Trial Info
    trialRIPlog = false(size(swrWindows,1),1);
    trialRIPlogPRE = false(size(swrWindows,1),1);
    trialRIPlogTRIAL = false(size(swrWindows,1),1);
    trialRIPlogPOST = false(size(swrWindows,1),1);
    trialsWithRips = nan(size(sessionMatrix,1),1);
    trialSWRlat = cell(size(sessionMatrix,1),1);
    trialSWRlatPRE = cell(size(sessionMatrix,1),1);
    trialSWRlatPOST = cell(size(sessionMatrix,1),1);
    for trl = 1:size(sessionMatrix,1)
        preTrlSWRlat = swrWindows(:,1) - sessionMatrix(trl,1);
        if trl == 1 || sessionMatrix(trl,4)==1
            preTrialSWRlog = (preTrlSWRlat<0) & (preTrlSWRlat>-1000);
        else
            preTrialSWRlog = (preTrlSWRlat<0) & (preTrlSWRlat>(sessionMatrix(trl-1,2)-sessionMatrix(trl,1))/2);
        end
        trialSWRlatPRE{trl} = preTrlSWRlat(preTrialSWRlog);
        trialRIPlogPRE(preTrialSWRlog) = true;
        
        trlSWRlog = (swrWindows(:,1) >= sessionMatrix(trl,1)) & (swrWindows(:,1) <= sessionMatrix(trl,2));
        trialSWRlat{trl} = swrWindows(trlSWRlog,1)-sessionMatrix(trl,1);
        trialsWithRips(trl) = sum(trlSWRlog) ~= 0;
        trialRIPlogTRIAL(trlSWRlog) = true;
        
        postTrlSWRlat = swrWindows(:,1) - sessionMatrix(trl,2);
        if trl == size(sessionMatrix,1) || sessionMatrix(trl,3)==0 || sessionMatrix(trl,4)==max(sessionMatrix(:,4))
            postTrialSWRlog = (postTrlSWRlat>0) & (postTrlSWRlat<1000);
        else
            postTrialSWRlog = (postTrlSWRlat>0) & (postTrlSWRlat<(sessionMatrix(trl+1,1)-sessionMatrix(trl,2))/2);
        end
        trialSWRlatPOST{trl} = postTrlSWRlat(postTrialSWRlog);
        trialRIPlogPOST(postTrialSWRlog) = true;
        
        trialRIPlog(preTrialSWRlog | trlSWRlog | postTrialSWRlog) = true;
    end
    
    % Non-Trial SWRs
    figure; 
    spDur = subplot(2,2,1);
    swrDurs = swrWindows(:,2)-swrWindows(:,1);
    histogram(swrDurs(~trialRIPlog), 0:10:max(swrDurs)+10, 'orientation', 'horizontal');
    set(spDur, 'XDir', 'reverse', 'YAxisLocation', 'right');
    title('SWR Duration');
    spPop = subplot(2,2,4);
    histogram(swrActv(~trialRIPlog), 0:0.01:max(swrActv)+0.1);
    title('Proportion Active Neurons');
    spDurPopCorr = subplot(2,2,2);
    corrScatPlot(swrActv(~trialRIPlog)', swrDurs(~trialRIPlog), '% Active', 'Duration', [], [])
    title('Non-Trial SWRs');
    linkaxes([spDur spDurPopCorr], 'y');
    linkaxes([spPop spDurPopCorr], 'x');
    supIRI = subplot(2,2,3);
    swrIRI = swrWindows(2:end,2)-swrWindows(1:end-1,1);
    histogram(swrIRI);
    title('Inter-SWR-Interval');
    annotation(gcf,'textbox', [0 0.95 1 0.05],'String', sprintf('%s Non-Trial SWRs', aniInfo(a).ID),...
        'FontWeight', 'Bold', 'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
    
    % Trial Period SWRs
    figure('toolbar', 'none');
    spDurPRE = subplot(6,2,1);
    histogram(swrDurs(trialRIPlogPRE), 0:50:max(swrDurs)+50, 'orientation', 'horizontal');
    set(spDurPRE, 'XDir', 'reverse', 'YAxisLocation', 'right');
    title('SWR Duration');
    spDurPopCorrPRE = subplot(6,2,2);
    corrScatPlot(swrActv(trialRIPlogPRE)', swrDurs(trialRIPlogPRE), '% Active', 'Duration', [], []);
    title('Pre Trial SWRs');
    spPRElat = subplot(6,2,3);
    histogram(cell2mat(trialSWRlatPRE), -4000:200:0);
    title('Pre Trial Latency');
    spPopPRE = subplot(6,2,4);
    histogram(swrActv(trialRIPlogPRE), 0:0.05:max(swrActv)+0.05);
    title('Proportion Active Neurons');    
    
    spDurTRL = subplot(6,2,5);
    histogram(swrDurs(trialRIPlogTRIAL), 0:50:max(swrDurs)+50, 'orientation', 'horizontal');
    set(spDurTRL, 'XDir', 'reverse', 'YAxisLocation', 'right');
    title('SWR Duration');
    spDurPopCorrTRL = subplot(6,2,6);
    corrScatPlot(swrActv(trialRIPlogTRIAL)', swrDurs(trialRIPlogTRIAL), '% Active', 'Duration', [], []);
    title('Trial Period SWRs');
    spTRLlat = subplot(6,2,7);
    histogram(cell2mat(trialSWRlat), 0:200:max(cell2mat(trialSWRlat)));
    title('Trial Latency');
    spPopTRL = subplot(6,2,8);
    histogram(swrActv(trialRIPlogTRIAL), 0:0.05:max(swrActv)+0.05);
    title('Proportion Active Neurons');    
    
    spDurPOST = subplot(6,2,9);
    histogram(swrDurs(trialRIPlogPOST), 0:50:max(swrDurs)+50, 'orientation', 'horizontal');
    set(spDurPOST, 'XDir', 'reverse', 'YAxisLocation', 'right');
    title('SWR Duration');
    spDurPopCorrPOST = subplot(6,2,10);
    corrScatPlot(swrActv(trialRIPlogPOST)', swrDurs(trialRIPlogPOST), '% Active', 'Duration', [], []);
    title('Post Trial SWRs');
    spPOSTlat = subplot(6,2,11);
    histogram(cell2mat(trialSWRlatPOST), 0:200:4000);
    title('Post Trial Latency');
    spPopPOST = subplot(6,2,12);
    histogram(swrActv(trialRIPlogPOST), 0:0.05:max(swrActv)+0.05);
    title('Proportion Active Neurons');    
       
    linkaxes([spDurPRE spDurTRL spDurPOST spDurPopCorrPRE spDurPopCorrTRL spDurPopCorrPOST], 'y');
    linkaxes([spPopPRE spPopTRL spPopPOST spDurPopCorrPRE spDurPopCorrTRL spDurPopCorrPOST], 'x');
    
    linkaxes([spPRElat, spTRLlat, spPOSTlat], 'y');  
    annotation(gcf,'textbox', [0 0.95 1 0.05],'String', sprintf('%s Trial SWRs', aniInfo(a).ID),...
        'FontWeight', 'Bold', 'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
    
    %% Calculate Trial Templates to evaluate behavioral relations of SWR events
    perfLog = sessionMatrix(:,3)==1;
    isLog = sessionMatrix(:,4)==sessionMatrix(:,5);
    
    templates = cell(1,max(sessionMatrix(:,5)));
    for o = 1:max(sessionMatrix(:,5))
        trlLog = sessionMatrix(:,5) == o & perfLog & isLog;
        tempSsnMtx = sessionMatrix(trlLog,:);
        tempTemplate = nan(240,size(ensemble.ensembleMatrix,2)-1, size(tempSsnMtx,1));
        for trl = 1:size(tempSsnMtx,1)
            tempTrial = ensemble.ensembleMatrix(tempSsnMtx(trl,1):tempSsnMtx(trl,1)+1199,2:end);
            for uni = 1:size(tempTrial,2)
                tempConv = conv(tempTrial(:,uni), ones(1,10)./(10/1000),'same');
                tempTemplate(:,uni,trl) = downsample(tempConv,5);
            end            
        end
        templates{o} = mean(tempTemplate,3);
    end
    catTemplate = cell2mat(templates');
    timeNdx = downsample(repmat(1:1200, [1,5]), 5);
    odrNdx = downsample([ones(1,1200), ones(1,1200)+1, ones(1,1200)+2, ones(1,1200)+3, ones(1,1200)+4], 5);
    
    %% All SWRs
    trlRelRips = swrWindows;
    ripDecodeODR = cell(size(trlRelRips,1),1);
    ripDecodeTIME = cell(size(trlRelRips,1),1);
    for rip = 1:size(trlRelRips,1)
        tempEnsemble = ensemble.ensembleMatrix(trlRelRips(rip,1):trlRelRips(rip,2),2:end);
        tempEnsConv = nan(size(tempEnsemble));
        for uni = 1:size(tempEnsemble,2)
            tempEnsConv(:,uni) = conv(tempEnsemble(:,uni), ones(1,10)./(10/1000), 'same');
        end
        post = CalcStaticBayesPost(catTemplate, downsample(tempEnsConv,10), 10);                 % <---------------- SWR Observation Rate
        ripDecodeTIME{rip} = DecodeBayesPost(post, timeNdx);
        ripDecodeODR{rip} = DecodeBayesPost(post, odrNdx);
    end
    
    lagBins = 0:1:100;
    odrLagDecodeMtx = repmat({zeros(max(sessionMatrix(:,4)), length(lagBins)-1)}, [1, max(sessionMatrix(:,4))]);
    odrLagRips = cell(length(ripDecodeODR), max(sessionMatrix(:,4)));
    tempHistBins = -1200:100:1200;
    timeLag = zeros(1,length(tempHistBins)-1);
    odrLagTrans = zeros(max(sessionMatrix(:,4)));
    for rip = 1:length(ripDecodeODR)
        tempRipODR = ripDecodeODR{rip};
        tempRipTIME = ripDecodeTIME{rip};
        tempRipNdx = 1:length(tempRipODR);
        tempODRlagDecodeMtx = odrLagDecodeMtx;
        for t = 1:length(tempRipODR)
            curOdr = tempRipODR(t);
            if ~isnan(curOdr)
                for o = 1:max(sessionMatrix(:,4))
                    odrTargLog = tempRipODR==o;
                    odrTargLag = tempRipNdx(odrTargLog) - t;
                    tempODRlagDecodeMtx{curOdr}(o,:) = tempODRlagDecodeMtx{curOdr}(o,:) + histcounts(odrTargLag, lagBins);
                end
            end
            if t<length(tempRipODR) && ~isnan(curOdr) && ~isnan(tempRipODR(t+1))
                odrLagTrans(curOdr,tempRipODR(t+1)) = odrLagTrans(curOdr,tempRipODR(t+1)) + 1;
                timeLag = timeLag + histcounts(tempRipTIME(t+1)-tempRipTIME(t), tempHistBins);
            end
        end
        odrLagRips(rip,:) = tempODRlagDecodeMtx;
    end
    odrLagNorm = nan(max(sessionMatrix(:,4)), length(lagBins)-1, max(sessionMatrix(:,4)));
    for o = 1:max(sessionMatrix(:,4))
        odrLagNorm(:,:,o) = mean(cell2mat(reshape(odrLagRips(:,o), [1,1,size(odrLagRips,1)])),3);
    end
    
    figure;
    subplot(5,1,1);
    aPlot = plot(odrLagNorm(:,:,1)');
    for o = 1:5
        aPlot(o).Color = odorColors(o,:);
    end
    title([aniInfo(a).ID ' All SWRs']);
    subplot(5,1,2);
    bPlot = plot(odrLagNorm(:,:,2)');
    for o = 1:5
        bPlot(o).Color = odorColors(o,:);
    end
    subplot(5,1,3);
    cPlot = plot(odrLagNorm(:,:,3)');
    for o = 1:5
        cPlot(o).Color = odorColors(o,:);
    end
    subplot(5,1,4);
    dPlot = plot(odrLagNorm(:,:,4)');
    for o = 1:5
        dPlot(o).Color = odorColors(o,:);
    end
    subplot(5,1,5);
    ePlot = plot(odrLagNorm(:,:,5)');
    for o = 1:5
        ePlot(o).Color = odorColors(o,:);
    end
    figure;
    subplot(2,2,1)
    imagesc(odrLagTrans);
    title(aniInfo(a).ID);
    subplot(2,2,2)
    bar(tempHistBins(1:end-1)+0.25,timeLag)
    set(gca, 'yscale', 'log')
    %% Pre Trial SWRs
    trlRelRips = swrWindows(trialRIPlogPRE,:);
%     trlRelRips = swrWindows;
    ripDecodeODR = cell(size(trlRelRips,1),1);
    ripDecodeTIME = cell(size(trlRelRips,1),1);
    for rip = 1:size(trlRelRips,1)
        tempEnsemble = ensemble.ensembleMatrix(trlRelRips(rip,1):trlRelRips(rip,2),2:end);
        tempEnsConv = nan(size(tempEnsemble));
        for uni = 1:size(tempEnsemble,2)
            tempEnsConv(:,uni) = conv(tempEnsemble(:,uni), ones(1,10)./(10/1000), 'same');
        end
        post = CalcStaticBayesPost(catTemplate, downsample(tempEnsConv,10), 10);                 % <---------------- SWR Observation Rate
        ripDecodeTIME{rip} = DecodeBayesPost(post, timeNdx);
        ripDecodeODR{rip} = DecodeBayesPost(post, odrNdx);
    end
    
    lagBins = 0:1:100;
    odrLagDecodeMtx = repmat({zeros(max(sessionMatrix(:,4)), length(lagBins)-1)}, [1, max(sessionMatrix(:,4))]);
    odrLagRips = cell(length(ripDecodeODR), max(sessionMatrix(:,4)));
    tempHistBins = -1200:100:1200;
    timeLag = zeros(1,length(tempHistBins)-1);
    odrLagTrans = zeros(max(sessionMatrix(:,4)));
    for rip = 1:length(ripDecodeODR)
        tempRipODR = ripDecodeODR{rip};
        tempRipTIME = ripDecodeTIME{rip};
        tempRipNdx = 1:length(tempRipODR);
        tempODRlagDecodeMtx = odrLagDecodeMtx;
        for t = 1:length(tempRipODR)
            curOdr = tempRipODR(t);
            if ~isnan(curOdr)
                for o = 1:max(sessionMatrix(:,4))
                    odrTargLog = tempRipODR==o;
                    odrTargLag = tempRipNdx(odrTargLog) - t;
                    tempODRlagDecodeMtx{curOdr}(o,:) = tempODRlagDecodeMtx{curOdr}(o,:) + histcounts(odrTargLag, lagBins);
                end
            end
            if t<length(tempRipODR) && ~isnan(curOdr) && ~isnan(tempRipODR(t+1))
                odrLagTrans(curOdr,tempRipODR(t+1)) = odrLagTrans(curOdr,tempRipODR(t+1)) + 1;
                timeLag = timeLag + histcounts(tempRipTIME(t+1)-tempRipTIME(t), tempHistBins);
            end
        end
        odrLagRips(rip,:) = tempODRlagDecodeMtx;
    end
    odrLagNorm = nan(max(sessionMatrix(:,4)), length(lagBins)-1, max(sessionMatrix(:,4)));
    for o = 1:max(sessionMatrix(:,4))
        odrLagNorm(:,:,o) = mean(cell2mat(reshape(odrLagRips(:,o), [1,1,size(odrLagRips,1)])),3);
    end
    
    figure;
    subplot(5,1,1);
    aPlot = plot(odrLagNorm(:,:,1)');
    for o = 1:5
        aPlot(o).Color = odorColors(o,:);
    end
    title([aniInfo(a).ID, ' Pre-Trial SWRs']);
    subplot(5,1,2);
    bPlot = plot(odrLagNorm(:,:,2)');
    for o = 1:5
        bPlot(o).Color = odorColors(o,:);
    end
    subplot(5,1,3);
    cPlot = plot(odrLagNorm(:,:,3)');
    for o = 1:5
        cPlot(o).Color = odorColors(o,:);
    end
    subplot(5,1,4);
    dPlot = plot(odrLagNorm(:,:,4)');
    for o = 1:5
        dPlot(o).Color = odorColors(o,:);
    end
    subplot(5,1,5);
    ePlot = plot(odrLagNorm(:,:,5)');
    for o = 1:5
        ePlot(o).Color = odorColors(o,:);
    end
    
    %% Post Trial SWRs
    trlRelRips = swrWindows(trialRIPlogPOST,:);
%     trlRelRips = swrWindows;
    ripDecodeODR = cell(size(trlRelRips,1),1);
    ripDecodeTIME = cell(size(trlRelRips,1),1);
    for rip = 1:size(trlRelRips,1)
        tempEnsemble = ensemble.ensembleMatrix(trlRelRips(rip,1):trlRelRips(rip,2),2:end);
        tempEnsConv = nan(size(tempEnsemble));
        for uni = 1:size(tempEnsemble,2)
            tempEnsConv(:,uni) = conv(tempEnsemble(:,uni), ones(1,10)./(10/1000), 'same');
        end
        post = CalcStaticBayesPost(catTemplate, downsample(tempEnsConv,10), 10);                 % <---------------- SWR Observation Rate
        ripDecodeTIME{rip} = DecodeBayesPost(post, timeNdx);
        ripDecodeODR{rip} = DecodeBayesPost(post, odrNdx);
    end
    
    lagBins = 0:1:100;
    odrLagDecodeMtx = repmat({zeros(max(sessionMatrix(:,4)), length(lagBins)-1)}, [1, max(sessionMatrix(:,4))]);
    odrLagRips = cell(length(ripDecodeODR), max(sessionMatrix(:,4)));
    tempHistBins = -1200:100:1200;
    timeLag = zeros(1,length(tempHistBins)-1);
    odrLagTrans = zeros(max(sessionMatrix(:,4)));
    for rip = 1:length(ripDecodeODR)
        tempRipODR = ripDecodeODR{rip};
        tempRipTIME = ripDecodeTIME{rip};
        tempRipNdx = 1:length(tempRipODR);
        tempODRlagDecodeMtx = odrLagDecodeMtx;
        for t = 1:length(tempRipODR)
            curOdr = tempRipODR(t);
            if ~isnan(curOdr)
                for o = 1:max(sessionMatrix(:,4))
                    odrTargLog = tempRipODR==o;
                    odrTargLag = tempRipNdx(odrTargLog) - t;
                    tempODRlagDecodeMtx{curOdr}(o,:) = tempODRlagDecodeMtx{curOdr}(o,:) + histcounts(odrTargLag, lagBins);
                end
            end
            if t<length(tempRipODR) && ~isnan(curOdr) && ~isnan(tempRipODR(t+1))
                odrLagTrans(curOdr,tempRipODR(t+1)) = odrLagTrans(curOdr,tempRipODR(t+1)) + 1;
                timeLag = timeLag + histcounts(tempRipTIME(t+1)-tempRipTIME(t), tempHistBins);
            end
        end
        odrLagRips(rip,:) = tempODRlagDecodeMtx;
    end
    odrLagNorm = nan(max(sessionMatrix(:,4)), length(lagBins)-1, max(sessionMatrix(:,4)));
    for o = 1:max(sessionMatrix(:,4))
        odrLagNorm(:,:,o) = mean(cell2mat(reshape(odrLagRips(:,o), [1,1,size(odrLagRips,1)])),3);
    end
    
    figure;
    subplot(5,1,1);
    aPlot = plot(odrLagNorm(:,:,1)');
    for o = 1:5
        aPlot(o).Color = odorColors(o,:);
    end
    title([aniInfo(a).ID ' Post-Trial SWRs']);
    subplot(5,1,2);
    bPlot = plot(odrLagNorm(:,:,2)');
    for o = 1:5
        bPlot(o).Color = odorColors(o,:);
    end
    subplot(5,1,3);
    cPlot = plot(odrLagNorm(:,:,3)');
    for o = 1:5
        cPlot(o).Color = odorColors(o,:);
    end
    subplot(5,1,4);
    dPlot = plot(odrLagNorm(:,:,4)');
    for o = 1:5
        dPlot(o).Color = odorColors(o,:);
    end
    subplot(5,1,5);
    ePlot = plot(odrLagNorm(:,:,5)');
    for o = 1:5
        ePlot(o).Color = odorColors(o,:);
    end
    %% Non Trial SWRs
    trlRelRips = swrWindows(~trialRIPlog,:);
%     trlRelRips = swrWindows;
    ripDecodeODR = cell(size(trlRelRips,1),1);
    ripDecodeTIME = cell(size(trlRelRips,1),1);
    for rip = 1:size(trlRelRips,1)
        tempEnsemble = ensemble.ensembleMatrix(trlRelRips(rip,1):trlRelRips(rip,2),2:end);
        tempEnsConv = nan(size(tempEnsemble));
        for uni = 1:size(tempEnsemble,2)
            tempEnsConv(:,uni) = conv(tempEnsemble(:,uni), ones(1,10)./(10/1000), 'same');
        end
        post = CalcStaticBayesPost(catTemplate, downsample(tempEnsConv,10), 10);                 % <---------------- SWR Observation Rate
        ripDecodeTIME{rip} = DecodeBayesPost(post, timeNdx);
        ripDecodeODR{rip} = DecodeBayesPost(post, odrNdx);
    end
    
    lagBins = 0:1:100;
    odrLagDecodeMtx = repmat({zeros(max(sessionMatrix(:,4)), length(lagBins)-1)}, [1, max(sessionMatrix(:,4))]);
    odrLagRips = cell(length(ripDecodeODR), max(sessionMatrix(:,4)));
    tempHistBins = -1200:100:1200;
    timeLag = zeros(1,length(tempHistBins)-1);
    odrLagTrans = zeros(max(sessionMatrix(:,4)));
    for rip = 1:length(ripDecodeODR)
        tempRipODR = ripDecodeODR{rip};
        tempRipTIME = ripDecodeTIME{rip};
        tempRipNdx = 1:length(tempRipODR);
        tempODRlagDecodeMtx = odrLagDecodeMtx;
        for t = 1:length(tempRipODR)
            curOdr = tempRipODR(t);
            if ~isnan(curOdr)
                for o = 1:max(sessionMatrix(:,4))
                    odrTargLog = tempRipODR==o;
                    odrTargLag = tempRipNdx(odrTargLog) - t;
                    tempODRlagDecodeMtx{curOdr}(o,:) = tempODRlagDecodeMtx{curOdr}(o,:) + histcounts(odrTargLag, lagBins);
                end
            end
            if t<length(tempRipODR) && ~isnan(curOdr) && ~isnan(tempRipODR(t+1))
                odrLagTrans(curOdr,tempRipODR(t+1)) = odrLagTrans(curOdr,tempRipODR(t+1)) + 1;
                timeLag = timeLag + histcounts(tempRipTIME(t+1)-tempRipTIME(t), tempHistBins);
            end
        end
        odrLagRips(rip,:) = tempODRlagDecodeMtx;
    end
    odrLagNorm = nan(max(sessionMatrix(:,4)), length(lagBins)-1, max(sessionMatrix(:,4)));
    for o = 1:max(sessionMatrix(:,4))
        odrLagNorm(:,:,o) = mean(cell2mat(reshape(odrLagRips(:,o), [1,1,size(odrLagRips,1)])),3);
    end
    
    figure;
    subplot(5,1,1);
    aPlot = plot(odrLagNorm(:,:,1)');
    for o = 1:5
        aPlot(o).Color = odorColors(o,:);
    end
    title([aniInfo(a).ID ' Non-Trial SWRs']);
    subplot(5,1,2);
    bPlot = plot(odrLagNorm(:,:,2)');
    for o = 1:5
        bPlot(o).Color = odorColors(o,:);
    end
    subplot(5,1,3);
    cPlot = plot(odrLagNorm(:,:,3)');
    for o = 1:5
        cPlot(o).Color = odorColors(o,:);
    end
    subplot(5,1,4);
    dPlot = plot(odrLagNorm(:,:,4)');
    for o = 1:5
        dPlot(o).Color = odorColors(o,:);
    end
    subplot(5,1,5);
    ePlot = plot(odrLagNorm(:,:,5)');
    for o = 1:5
        ePlot(o).Color = odorColors(o,:);
    end
    
    %% Plot the session SWR traces
    figure; 
    sp1 = subplot(2,1,1);
    plot(swTrace);
    hold on; 
    plot(ripTrace);
    title(aniInfo(a).ID);
      
    sp2 = subplot(2,1,2);
    plot(newSWtrace);
    hold on; 
    plot(newRIPtrace);
    
    linkaxes([sp1 sp2], 'xy');
    
    annotation(gcf,'textbox', [0 0.01 1 0.05],'String', ['SW File:' swFile],...
        'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
    annotation(gcf,'textbox', [0 0.05 1 0.05],'String', ['RIP File: ' ripFile],...
        'FontSize',10, 'edgecolor', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
end

%% %%%%%%%%%
%  Functions
%  %%%%%%%%%

%%
function [ripEpocs, lfpFilt, zPow] = RippleDetection(lfp, fs, powThresh, mergeThresh)
%% RippleDetection
% Code for detection of ripples on tetrode presumed to be in pyramidal
%   layer of CA1
% Inputs:
%   lfp: Nx1 input of the raw LFP 
%   powThresh: 1x2 input with threshold parameters for detection of peak
%       event. Initial value is event boundaries, second is peak required
%   mergeThresh: threshold for merging two sharpwave events if they happen
%       in close temporal proximity%
%%
if nargin == 0
    [file,path] = uigetfile('*.mat');
    load([path file]); %#ok<LOAD>
    ripCol = cellfun(@(a)~isempty(a), regexp(statMatrixColIDs, 'LFP_Ripple$'));
    lfp = statMatrix(:,ripCol); %#ok<USENS>
    clear statMatrix statMatrixColIDs file path
    powThresh = [0 4];
    mergeThresh = 15;
end
%%
lfpFilt = SimpleFilt(lfp, fs, 150, 250);
zPow = zscore(abs(hilbert(lfpFilt)));

%% Identify Potential Sharpwave Periods
abvThresh1 = zPow > powThresh(1);
epocStart = find(diff(abvThresh1)==1);
epocEnd = find(diff(abvThresh1)==-1);
epocWindows = [epocStart, nan(size(epocStart))];
for epoc = 1:size(epocStart,1)
    tempEpocEnd = epocEnd(find(epocEnd>epocStart(epoc), 1, 'first'));
    if sum(zPow(epocStart(epoc):tempEpocEnd)>=powThresh(2))>=1
        epocWindows(epoc,2) = tempEpocEnd;
    else
        epocWindows(epoc,1) = nan;
    end
end
epocWindows(isnan(epocWindows(:,1)),:) = [];
% Merge short latency sharpwaves
for epoc = 2:size(epocWindows,1)
    if ~isnan(epocWindows(epoc,1)) && ~isnan(epocWindows(epoc-1,1))
        if epocWindows(epoc,1) - epocWindows(epoc-1,2) <= mergeThresh
            epocWindows(epoc,1) = epocWindows(epoc-1,1);
            epocWindows(epoc-1,:) = nan;
        end
    end
end
epocWindows(isnan(epocWindows(:,1)),:) = [];
%%
ripEpocs = epocWindows;
end

function [swEpocs, hpf, zPow] = SharpwaveDetection(lfp, fs, powThresh, mergeThresh)
%% SharpwaveDetection
% Code for detection of sharpwaves on tetrodes presumed to be in the 
%   stratum radiatum of CA1
% Inputs:
%   lfp: Nx1 input of the raw LFP 
%   powThresh: 1x2 input with threshold parameters for detection of peak
%       event. Initial value is event boundaries, second is peak required
%   mergeThresh: threshold for merging two sharpwave events if they happen
%       in close temporal proximity
%%
if nargin == 0
    [file,path] = uigetfile('*.mat');
    load([path file]); %#ok<LOAD>
    lfp = statMatrix(:,2); %#ok<USENS>
    fs = 1/mode(diff(statMatrix(:,1)));
    clear statMatrix statMatrixColIDs file path
    powThresh = [0 2];
    mergeThresh = 15;
end
%% Filter & Z-Score Envelope
hpf = highpass(lfp, 4, fs);
zPow = zscore(abs(hilbert(hpf)));

%% Identify Potential Sharpwave Periods
abvThresh1 = zPow > powThresh(1);
epocStart = find(diff(abvThresh1)==1);
epocEnd = find(diff(abvThresh1)==-1);
epocWindows = [epocStart, nan(size(epocStart))];
for epoc = 1:size(epocStart,1)
    tempEpocEnd = epocEnd(find(epocEnd>epocStart(epoc), 1, 'first'));
    if sum(zPow(epocStart(epoc):tempEpocEnd)>=powThresh(2))>=1
        epocWindows(epoc,2) = tempEpocEnd;
    else
        epocWindows(epoc,1) = nan;
    end
end
epocWindows(isnan(epocWindows(:,1)),:) = [];
% Merge short latency sharpwaves
for epoc = 2:size(epocWindows,1)
    if ~isnan(epocWindows(epoc,1)) && ~isnan(epocWindows(epoc-1,1))
        if epocWindows(epoc,1) - epocWindows(epoc-1,2) <= mergeThresh
            epocWindows(epoc,1) = epocWindows(epoc-1,1);
            epocWindows(epoc-1,:) = nan;
        end
    end
end
epocWindows(isnan(epocWindows(:,1)),:) = [];
%%
swEpocs = epocWindows;
%%
end

%% Simple Filtering
function lfpFilt = SimpleFilt(trace, samp, low, high)
Wn_FRange = [low/(samp/2) high/(samp/2)]; % normalized by the nyquist frequency
[bFRange, aFRange] = butter(3,Wn_FRange);
lfpFilt = filtfilt(bFRange,aFRange,trace);
end

%%
function post = CalcStaticBayesPost(likely, obsv, binSize)
post = nan(size(obsv,1), size(likely,1), size(obsv,3));
for trl = 1:size(obsv,3)
    for t = 1:size(obsv,1)
        p = nan(size(likely));
        curPopVect = obsv(t,:,trl)*(binSize/1000);
        curPopFact = factorial(curPopVect);
        for u = 1:size(likely,2)
            curAvgUniFR = likely(:,u);
            p(:,u) = (((binSize/1000).*curAvgUniFR).^curPopVect(u))./curPopFact(u);
        end        
        pp = prod(p,2);
        ee = exp(-((binSize/1000)*sum(likely,2)));
        tempPost = pp.*ee;
        post(t,:,trl) = tempPost./sum(tempPost);
    end
end
end

%% 
function decode = DecodeBayesPost(post, id)
% Assumes post is in the structure of ObservTime X LikelyTime X Trial
decode = nan(size(post,1),size(post,3));
for o = 1:size(post,3)
    for d = 1:size(post,1)
        if ~isnan(post(d,1,o))
            decode(d,o) = id(find(post(d,:,o)==max(post(d,:,o)),1,'first'));
        end
    end
end
end