ticMain = tic;
%% Variables
terminalIndex = 2472289;
windowSize = 1000;

thetaWndw = [4 12];
lowBetaWndw = [13 19];
betaWndw = [20 40];
lowGammaWndw = [41 59];
highGammaWndw = [61 80];
rippleWndw = [150 250];

%% Load Directory and shiet
oefileDir = uigetdir;
oeFileDirFiles = dir(oefileDir);
oeFiles = {oeFileDirFiles.name};
fpgaFiles = oeFiles(logical(cellfun(@(a)~isempty(a), strfind(oeFiles, '100'))));
lfpFiles = fpgaFiles(logical(cellfun(@(a)~isempty(a), strfind(fpgaFiles, 'CH'))));
ssnDataFile = oeFiles{logical(cellfun(@(a)~isempty(a), strfind(oeFiles, '.mat')))};

oeSessionInfo = get_session_info(oefileDir);
sampleRate = oeSessionInfo.sampleRate;

%% Extract OpenEphys Event Timestamps
[eventIndexVects,onTS,offTS,eventTS] = BehavTS(oefileDir,100,0);

% Create Poke Events Column 
% NOTE: This is an unedited vector, i.e. it has ALL the poke events in it,
% even ones that were originally excluded as "short pokes" during data
% collection by the matlab code running the session. DO NOT USE THIS CODE
% WITHOUT ACCOUNTING FOR THAT!!! (i.e. don't just use this vector to align
% things to, make sure you're pulling the trial relevant poke events).
pokeInTSs = onTS.OnPeakTS_5;
pokeInNdxs = nan(length(pokeInTSs),1);
for pI = 1:length(pokeInTSs)
    pokeInNdxs(pI) = find(pokeInTSs(pI)>=eventTS, 1, 'last');
end
pokeOutTSs = offTS.OffPeakTS_5;
pokeOutNdxs = nan(length(pokeOutTSs),1);
for pO = 1:length(pokeOutTSs)
    pokeOutNdxs(pO) = find(pokeOutTSs(pO)>=eventTS, 1, 'last');
end    

% Pull Out Odor Delivery Order & Odor Delivery Indices
% Odor A
odorAonTSs = onTS.OnPeakTS_1;
odorAoffTSs = offTS.OffPeakTS_1;
odorAndxs = nan(length(odorAonTSs),3);
for oA = 1:length(odorAndxs)
    odorAndxs(oA,1) = find(eventTS<=odorAonTSs(oA), 1, 'last');
    odorAndxs(oA,2) = 1;
    if oA>length(odorAoffTSs)
        odorAndxs(oA,3) = nan;
    else
        odorAndxs(oA,3) = find(eventTS<=odorAoffTSs(oA), 1, 'last');
    end
end

odorBonTSs = onTS.OnPeakTS_2;
odorBoffTSs = offTS.OffPeakTS_2;
odorBndxs = nan(length(odorBonTSs),3);
for oB = 1:length(odorBndxs)
    odorBndxs(oB,1) = find(eventTS<=odorBonTSs(oB), 1, 'last');
    odorBndxs(oB,2) = 2;
    if oB>length(odorBoffTSs)
        odorBndxs(oB,3) = nan;
    else
        odorBndxs(oB,3) = find(eventTS<=odorBoffTSs(oB), 1, 'last');
    end
end

odorConTSs = onTS.OnPeakTS_3;
odorCoffTSs = offTS.OffPeakTS_3;
odorCndxs = nan(length(odorConTSs),3);
for oC = 1:length(odorCndxs)
    odorCndxs(oC,1) = find(eventTS<=odorConTSs(oC), 1, 'last');
    odorCndxs(oC,2) = 3;
    if oC>length(odorCoffTSs)
        odorCndxs(oC,3) = nan;
    else
        odorCndxs(oC,3) = find(eventTS<=odorCoffTSs(oC), 1, 'last');
    end
end

odorDonTSs = onTS.OnPeakTS_4;
odorDoffTSs = offTS.OffPeakTS_4;
odorDndxs = nan(length(odorDonTSs),3);
for oD = 1:length(odorDndxs)
    odorDndxs(oD,1) = find(eventTS<=odorDonTSs(oD), 1, 'last');
    odorDndxs(oD,2) = 4;
    if oD>length(odorDoffTSs)
        odorDndxs(oD,3) = nan;
    else
        odorDndxs(oD,3) = find(eventTS<=odorDoffTSs(oD), 1, 'last');
    end
end

oeOdorIDnNDXs = sortrows([odorAndxs; odorBndxs; odorCndxs; odorDndxs]);

% Check for duplicate odor entries... this happens when someone forgets to
% stop recording and depressurizes the system after a behavior session.
if length(unique(oeOdorIDnNDXs(:,1))) < size(oeOdorIDnNDXs,1)
    duplicatesToRemove = false(size(oeOdorIDnNDXs,1),1);
    for o = 1:size(oeOdorIDnNDXs,1)
        if sum(oeOdorIDnNDXs(:,1)==oeOdorIDnNDXs(o,1))>=2
            duplicatesToRemove(oeOdorIDnNDXs(:,1)==oeOdorIDnNDXs(o,1)) = true;
        end
    end
    oeOdorIDnNDXs(duplicatesToRemove,:) = [];
end            

%% Extract ssnData Variables
load([oefileDir '\' ssnDataFile]);

% You'll notice I step through the ssnData structure trial by trial here to
% pull out poke in/out indices, and then again below in the next section. I
% could forgo the step here and just do it down below, but I'd rather do it
% separately here, before I start loading files and doing all the
% filtering/power extraction etc because that takes significantly more time
% than doing this quick check before loading any of the LFP data.
pokeInTrlNdx = nan(length(ssnData),1);
pokeOutTrlNdx = nan(length(ssnData),1);
for trl = 1:length(ssnData)
    curTrlOdor = ssnData(trl).Odor;
    if curTrlOdor == oeOdorIDnNDXs(trl,2)
        curTrlInSeqLog = ssnData(trl).TranspositionDistance==0;
        curTrlPerfLog = ssnData(trl).Performance==1;
        
        pokeInTrlNdx(trl) = pokeInNdxs(find(pokeInNdxs<=oeOdorIDnNDXs(trl,1),1,'last'));
        if (curTrlInSeqLog && curTrlPerfLog) || (~curTrlInSeqLog && ~curTrlPerfLog)
            pokeOutTrlNdx(trl) = pokeOutNdxs(find(pokeOutNdxs>=oeOdorIDnNDXs(trl,3),1,'first'));
        elseif (curTrlInSeqLog && ~curTrlPerfLog) || (~curTrlInSeqLog && curTrlPerfLog)
            pokeOutTrlNdx(trl) = pokeOutNdxs(find(pokeOutNdxs<=oeOdorIDnNDXs(trl,3),1,'last'));
        end
    else
        error('Odor identity does not match for trial %i. OpenEphys = %i; ssnData = %i. Check files and correct data sources as necessary', trl, curTrlOdor, oeOdorIDnNDXs(trl,2));
    end
end

if sum(pokeInTrlNdx>(terminalIndex-(windowSize/2)))>=1
    lastGoodTrl = find(pokeInTrlNdx<(terminalIndex-(windowSize/2)),1,'last');
    ssnData = ssnData(1:lastGoodTrl);
end

%% Compile the data
thetaTrlPower = cell(length(ssnData), length(lfpFiles),2);                  thetaTrlPhase = cell(length(ssnData), length(lfpFiles),2);
lowBetaTrlPower = cell(length(ssnData), length(lfpFiles),2);                lowBetaTrlPhase = cell(length(ssnData), length(lfpFiles),2);
betaTrlPower = cell(length(ssnData), length(lfpFiles),2);                   betaTrlPhase = cell(length(ssnData), length(lfpFiles),2);
lowGammaTrlPower = cell(length(ssnData), length(lfpFiles),2);               lowGammaTrlPhase = cell(length(ssnData), length(lfpFiles),2);
highGammaTrlPower = cell(length(ssnData), length(lfpFiles),2);              highGammaTrlPhase = cell(length(ssnData), length(lfpFiles),2);
rippleTrlPower = cell(length(ssnData), length(lfpFiles),2);                 rippleTrlPhase = cell(length(ssnData), length(lfpFiles),2);

tic
for ch = 1:length(lfpFiles)    
    [LFP_Raw,~,~] = load_open_ephys_data_faster([oefileDir '\' lfpFiles{ch}]);
    LFP_Raw = LFP_Raw(1:terminalIndex);
    eventTS = eventTS(1:terminalIndex);
    
    [thetaHilb, thetaFilt] = PhaseFreqDetectAbbr(LFP_Raw, eventTS, thetaWndw(1), thetaWndw(2));    
    [thetaRMS,~] = mvrms_mfile_simulink(eventTS, thetaFilt, 1.25*mean(thetaWndw), 0);
    thetaRMS = zscore(thetaRMS);
    
    [lowBetaHilb, lowBetaFilt] = PhaseFreqDetectAbbr(LFP_Raw, eventTS, lowBetaWndw(1), lowBetaWndw(2));
    [lowBetaRMS,~] = mvrms_mfile_simulink(eventTS, lowBetaFilt, 1.25*mean(lowBetaWndw), 0);
    lowBetaRMS = zscore(lowBetaRMS);
    
    [betaHilb, betaFilt] = PhaseFreqDetectAbbr(LFP_Raw, eventTS, betaWndw(1), betaWndw(2));
    [betaRMS,~] = mvrms_mfile_simulink(eventTS, betaFilt, 1.25*mean(betaWndw), 0);
    betaRMS = zscore(betaRMS);
    
    [lowGammaHilb, lowGammaFilt] = PhaseFreqDetectAbbr(LFP_Raw, eventTS, lowGammaWndw(1), lowGammaWndw(2));
    [lowGammaRMS,~] = mvrms_mfile_simulink(eventTS, lowGammaFilt, 1.25*mean(lowGammaWndw), 0);
    lowGammaRMS = zscore(lowGammaRMS);
    
    [highGammaHilb, highGammaFilt] = PhaseFreqDetectAbbr(LFP_Raw, eventTS, highGammaWndw(1), highGammaWndw(2));
    [highGammaRMS,~] = mvrms_mfile_simulink(eventTS, highGammaFilt, 1.25*mean(highGammaWndw), 0);
    highGammaRMS = zscore(highGammaRMS);
    
    [rippleHilb, rippleFilt] = PhaseFreqDetectAbbr(LFP_Raw, eventTS, rippleWndw(1), rippleWndw(2));
    [rippleRMS,~] = mvrms_mfile_simulink(eventTS, rippleFilt, 1.25*mean(rippleWndw), 0);
    rippleRMS = zscore(rippleRMS);
    
    for trl = 1:length(ssnData)
        curPokeInNdx = pokeInTrlNdx(trl);
        curPokeOutNdx = pokeOutTrlNdx(trl);
        % Theta
        thetaTrlPower{trl,ch,1} = thetaRMS(curPokeInNdx-(windowSize/2):curPokeInNdx+(windowSize/2));
        thetaTrlPower{trl,ch,2} = thetaRMS(curPokeOutNdx-(windowSize/2):curPokeOutNdx+(windowSize/2));
        thetaTrlPhase{trl,ch,1} = thetaHilb(curPokeInNdx-(windowSize/2):curPokeInNdx+(windowSize/2));
        thetaTrlPhase{trl,ch,2} = thetaHilb(curPokeOutNdx-(windowSize/2):curPokeOutNdx+(windowSize/2));
        % LowBeta
        lowBetaTrlPower{trl,ch,1} = lowBetaRMS(curPokeInNdx-(windowSize/2):curPokeInNdx+(windowSize/2));
        lowBetaTrlPower{trl,ch,2} = lowBetaRMS(curPokeOutNdx-(windowSize/2):curPokeOutNdx+(windowSize/2));
        lowBetaTrlPhase{trl,ch,1} = lowBetaHilb(curPokeInNdx-(windowSize/2):curPokeInNdx+(windowSize/2));
        lowBetaTrlPhase{trl,ch,2} = lowBetaHilb(curPokeOutNdx-(windowSize/2):curPokeOutNdx+(windowSize/2));
        % Beta
        betaTrlPower{trl,ch,1} = betaRMS(curPokeInNdx-(windowSize/2):curPokeInNdx+(windowSize/2));
        betaTrlPower{trl,ch,2} = betaRMS(curPokeOutNdx-(windowSize/2):curPokeOutNdx+(windowSize/2));
        betaTrlPhase{trl,ch,1} = betaHilb(curPokeInNdx-(windowSize/2):curPokeInNdx+(windowSize/2));
        betaTrlPhase{trl,ch,2} = betaHilb(curPokeOutNdx-(windowSize/2):curPokeOutNdx+(windowSize/2));
        % LowGamma
        lowGammaTrlPower{trl,ch,1} = lowGammaRMS(curPokeInNdx-(windowSize/2):curPokeInNdx+(windowSize/2));
        lowGammaTrlPower{trl,ch,2} = lowGammaRMS(curPokeOutNdx-(windowSize/2):curPokeOutNdx+(windowSize/2));
        lowGammaTrlPhase{trl,ch,1} = lowGammaHilb(curPokeInNdx-(windowSize/2):curPokeInNdx+(windowSize/2));
        lowGammaTrlPhase{trl,ch,2} = lowGammaHilb(curPokeOutNdx-(windowSize/2):curPokeOutNdx+(windowSize/2));
        % HighGamma        
        highGammaTrlPower{trl,ch,1} = highGammaRMS(curPokeInNdx-(windowSize/2):curPokeInNdx+(windowSize/2));
        highGammaTrlPower{trl,ch,2} = highGammaRMS(curPokeOutNdx-(windowSize/2):curPokeOutNdx+(windowSize/2));
        highGammaTrlPhase{trl,ch,1} = highGammaHilb(curPokeInNdx-(windowSize/2):curPokeInNdx+(windowSize/2));
        highGammaTrlPhase{trl,ch,2} = highGammaHilb(curPokeOutNdx-(windowSize/2):curPokeOutNdx+(windowSize/2));
        % Ripple
        rippleTrlPower{trl,ch,1} = rippleRMS(curPokeInNdx-(windowSize/2):curPokeInNdx+(windowSize/2));
        rippleTrlPower{trl,ch,2} = rippleRMS(curPokeOutNdx-(windowSize/2):curPokeOutNdx+(windowSize/2));
        rippleTrlPhase{trl,ch,1} = rippleHilb(curPokeInNdx-(windowSize/2):curPokeInNdx+(windowSize/2));
        rippleTrlPhase{trl,ch,2} = rippleHilb(curPokeOutNdx-(windowSize/2):curPokeOutNdx+(windowSize/2));
    end
    fprintf('Compiled Channel %03i          ', ch);
    toc
end

%% Now run the analyses
thetaPowerTrlSimIn = nan(length(lfpFiles), length(lfpFiles), length(ssnData));          thetaPowerTrlSimOut = nan(length(lfpFiles), length(lfpFiles), length(ssnData));
thetaPhaseTrlSimIn = nan(length(lfpFiles), length(lfpFiles), length(ssnData));          thetaPhaseTrlSimOut = nan(length(lfpFiles), length(lfpFiles), length(ssnData));

lowBetaPowerTrlSimIn = nan(length(lfpFiles), length(lfpFiles), length(ssnData));        lowBetaPowerTrlSimOut = nan(length(lfpFiles), length(lfpFiles), length(ssnData));
lowBetaPhaseTrlSimIn = nan(length(lfpFiles), length(lfpFiles), length(ssnData));        lowBetaPhaseTrlSimOut = nan(length(lfpFiles), length(lfpFiles), length(ssnData));

betaPowerTrlSimIn = nan(length(lfpFiles), length(lfpFiles), length(ssnData));           betaPowerTrlSimOut = nan(length(lfpFiles), length(lfpFiles), length(ssnData));
betaPhaseTrlSimIn = nan(length(lfpFiles), length(lfpFiles), length(ssnData));           betaPhaseTrlSimOut = nan(length(lfpFiles), length(lfpFiles), length(ssnData));

lowGammaPowerTrlSimIn = nan(length(lfpFiles), length(lfpFiles), length(ssnData));       lowGammaPowerTrlSimOut = nan(length(lfpFiles), length(lfpFiles), length(ssnData));
lowGammaPhaseTrlSimIn = nan(length(lfpFiles), length(lfpFiles), length(ssnData));       lowGammaPhaseTrlSimOut = nan(length(lfpFiles), length(lfpFiles), length(ssnData));

highGammaPowerTrlSimIn = nan(length(lfpFiles), length(lfpFiles), length(ssnData));      highGammaPowerTrlSimOut = nan(length(lfpFiles), length(lfpFiles), length(ssnData));
highGammaPhaseTrlSimIn = nan(length(lfpFiles), length(lfpFiles), length(ssnData));      highGammaPhaseTrlSimOut = nan(length(lfpFiles), length(lfpFiles), length(ssnData));

ripplePowerTrlSimIn = nan(length(lfpFiles), length(lfpFiles), length(ssnData));         ripplePowerTrlSimOut = nan(length(lfpFiles), length(lfpFiles), length(ssnData));
ripplePhaseTrlSimIn = nan(length(lfpFiles), length(lfpFiles), length(ssnData));         ripplePhaseTrlSimOut = nan(length(lfpFiles), length(lfpFiles), length(ssnData));

tic
for trl = 1:length(ssnData)    
    thetaPowerTrlSimIn(:,:,trl) = corr(cell2mat(thetaTrlPower(trl,:,1)));
    thetaPowerTrlSimOut(:,:,trl) = corr(cell2mat(thetaTrlPower(trl,:,2)));
    
    lowBetaPowerTrlSimIn(:,:,trl) = corr(cell2mat(lowBetaTrlPower(trl,:,1)));
    lowBetaPowerTrlSimOut(:,:,trl) = corr(cell2mat(lowBetaTrlPower(trl,:,2)));
    
    betaPowerTrlSimIn(:,:,trl) = corr(cell2mat(betaTrlPower(trl,:,1)));
    betaPowerTrlSimOut(:,:,trl) = corr(cell2mat(betaTrlPower(trl,:,2)));
    
    lowGammaPowerTrlSimIn(:,:,trl) = corr(cell2mat(lowGammaTrlPower(trl,:,1)));
    lowGammaPowerTrlSimOut(:,:,trl) = corr(cell2mat(lowGammaTrlPower(trl,:,2)));
    
    highGammaPowerTrlSimIn(:,:,trl) = corr(cell2mat(highGammaTrlPower(trl,:,1)));
    highGammaPowerTrlSimOut(:,:,trl) = corr(cell2mat(highGammaTrlPower(trl,:,2)));
    
    ripplePowerTrlSimIn(:,:,trl) = corr(cell2mat(rippleTrlPower(trl,:,1)));
    ripplePowerTrlSimOut(:,:,trl) = corr(cell2mat(rippleTrlPower(trl,:,2)));
    
    for ch1 = 1:length(lfpFiles)
        for ch2 = 1:length(lfpFiles)
            [thetaPhaseTrlSimIn(ch1,ch2,trl), ~] = circ_corrcc(cell2mat(thetaTrlPhase(trl,ch1,1)), cell2mat(thetaTrlPhase(trl,ch2,1)));
            [thetaPhaseTrlSimOut(ch1,ch2,trl), ~] = circ_corrcc(cell2mat(thetaTrlPhase(trl,ch1,2)), cell2mat(thetaTrlPhase(trl,ch2,2)));
            
            [lowBetaPhaseTrlSimIn(ch1,ch2,trl), ~] = circ_corrcc(cell2mat(lowBetaTrlPhase(trl,ch1,1)), cell2mat(lowBetaTrlPhase(trl,ch2,1)));
            [lowBetaPhaseTrlSimOut(ch1,ch2,trl), ~] = circ_corrcc(cell2mat(lowBetaTrlPhase(trl,ch1,2)), cell2mat(lowBetaTrlPhase(trl,ch2,2)));
            
            [betaPhaseTrlSimIn(ch1,ch2,trl), ~] = circ_corrcc(cell2mat(betaTrlPhase(trl,ch1,1)), cell2mat(betaTrlPhase(trl,ch2,1)));
            [betaPhaseTrlSimOut(ch1,ch2,trl), ~] = circ_corrcc(cell2mat(betaTrlPhase(trl,ch1,2)), cell2mat(betaTrlPhase(trl,ch2,2)));
            
            [lowGammaPhaseTrlSimIn(ch1,ch2,trl), ~] = circ_corrcc(cell2mat(lowGammaTrlPhase(trl,ch1,1)), cell2mat(lowGammaTrlPhase(trl,ch2,1)));
            [lowGammaPhaseTrlSimOut(ch1,ch2,trl), ~] = circ_corrcc(cell2mat(lowGammaTrlPhase(trl,ch1,2)), cell2mat(lowGammaTrlPhase(trl,ch2,2)));
            
            [highGammaPhaseTrlSimIn(ch1,ch2,trl), ~] = circ_corrcc(cell2mat(highGammaTrlPhase(trl,ch1,1)), cell2mat(highGammaTrlPhase(trl,ch2,1)));
            [highGammaPhaseTrlSimOut(ch1,ch2,trl), ~] = circ_corrcc(cell2mat(highGammaTrlPhase(trl,ch1,2)), cell2mat(highGammaTrlPhase(trl,ch2,2)));
            
            [ripplePhaseTrlSimIn(ch1,ch2,trl), ~] = circ_corrcc(cell2mat(rippleTrlPhase(trl,ch1,1)), cell2mat(rippleTrlPhase(trl,ch2,1)));
            [ripplePhaseTrlSimOut(ch1,ch2,trl), ~] = circ_corrcc(cell2mat(rippleTrlPhase(trl,ch1,2)), cell2mat(rippleTrlPhase(trl,ch2,2)));
        end
    end
    fprintf('Analyzed Trial %02i          ', trl);
    toc
end
    
%% Separate Out Trial Types
inSeqCorrectLog = [ssnData.TranspositionDistance]==0 & [ssnData.Performance]==1;
inSeqIncorrectLog = [ssnData.TranspositionDistance]==0 & [ssnData.Performance]==0;
outSeqCorrectLog = ~([ssnData.TranspositionDistance]==0) & [ssnData.Performance]==1;
outSeqIncorrectLog = ~([ssnData.TranspositionDistance]==0) & [ssnData.Performance]==0;

% Plot Theta Power Similarity Matrices
PlotMeanPowerTrlSimilarityMatrix(thetaPowerTrlSimIn(:,:,inSeqCorrectLog), thetaPowerTrlSimOut(:,:,inSeqCorrectLog),...
    thetaPowerTrlSimIn(:,:,inSeqIncorrectLog), thetaPowerTrlSimOut(:,:,inSeqIncorrectLog),...
    thetaPowerTrlSimIn(:,:,outSeqCorrectLog), thetaPowerTrlSimOut(:,:,outSeqCorrectLog),...
    thetaPowerTrlSimIn(:,:,outSeqIncorrectLog), thetaPowerTrlSimOut(:,:,outSeqIncorrectLog), 'Theta-Power');
PlotMeanPowerTrlSimilarityMatrix(thetaPhaseTrlSimIn(:,:,inSeqCorrectLog), thetaPhaseTrlSimOut(:,:,inSeqCorrectLog),...
    thetaPhaseTrlSimIn(:,:,inSeqIncorrectLog), thetaPhaseTrlSimOut(:,:,inSeqIncorrectLog),...
    thetaPhaseTrlSimIn(:,:,outSeqCorrectLog), thetaPhaseTrlSimOut(:,:,outSeqCorrectLog),...
    thetaPhaseTrlSimIn(:,:,outSeqIncorrectLog), thetaPhaseTrlSimOut(:,:,outSeqIncorrectLog), 'Theta-Phase');

% Plot Low Beta Power Similarity Matrices
PlotMeanPowerTrlSimilarityMatrix(lowBetaPowerTrlSimIn(:,:,inSeqCorrectLog), lowBetaPowerTrlSimOut(:,:,inSeqCorrectLog),...
    lowBetaPowerTrlSimIn(:,:,inSeqIncorrectLog), lowBetaPowerTrlSimOut(:,:,inSeqIncorrectLog),...
    lowBetaPowerTrlSimIn(:,:,outSeqCorrectLog), lowBetaPowerTrlSimOut(:,:,outSeqCorrectLog),...
    lowBetaPowerTrlSimIn(:,:,outSeqIncorrectLog), lowBetaPowerTrlSimOut(:,:,outSeqIncorrectLog), 'LowBeta-Power');
PlotMeanPowerTrlSimilarityMatrix(lowBetaPhaseTrlSimIn(:,:,inSeqCorrectLog), lowBetaPhaseTrlSimOut(:,:,inSeqCorrectLog),...
    lowBetaPhaseTrlSimIn(:,:,inSeqIncorrectLog), lowBetaPhaseTrlSimOut(:,:,inSeqIncorrectLog),...
    lowBetaPhaseTrlSimIn(:,:,outSeqCorrectLog), lowBetaPhaseTrlSimOut(:,:,outSeqCorrectLog),...
    lowBetaPhaseTrlSimIn(:,:,outSeqIncorrectLog), lowBetaPhaseTrlSimOut(:,:,outSeqIncorrectLog), 'LowBeta-Phase');

% Plot Beta Power Similarity Matrices
PlotMeanPowerTrlSimilarityMatrix(betaPowerTrlSimIn(:,:,inSeqCorrectLog), betaPowerTrlSimOut(:,:,inSeqCorrectLog),...
    betaPowerTrlSimIn(:,:,inSeqIncorrectLog), betaPowerTrlSimOut(:,:,inSeqIncorrectLog),...
    betaPowerTrlSimIn(:,:,outSeqCorrectLog), betaPowerTrlSimOut(:,:,outSeqCorrectLog),...
    betaPowerTrlSimIn(:,:,outSeqIncorrectLog), betaPowerTrlSimOut(:,:,outSeqIncorrectLog), 'Beta-Power');
PlotMeanPowerTrlSimilarityMatrix(betaPhaseTrlSimIn(:,:,inSeqCorrectLog), betaPhaseTrlSimOut(:,:,inSeqCorrectLog),...
    betaPhaseTrlSimIn(:,:,inSeqIncorrectLog), betaPhaseTrlSimOut(:,:,inSeqIncorrectLog),...
    betaPhaseTrlSimIn(:,:,outSeqCorrectLog), betaPhaseTrlSimOut(:,:,outSeqCorrectLog),...
    betaPhaseTrlSimIn(:,:,outSeqIncorrectLog), betaPhaseTrlSimOut(:,:,outSeqIncorrectLog), 'Beta-Phase');

% Plot Low Gamma Power Similarity Matrices
PlotMeanPowerTrlSimilarityMatrix(lowGammaPowerTrlSimIn(:,:,inSeqCorrectLog), lowGammaPowerTrlSimOut(:,:,inSeqCorrectLog),...
    lowGammaPowerTrlSimIn(:,:,inSeqIncorrectLog), lowGammaPowerTrlSimOut(:,:,inSeqIncorrectLog),...
    lowGammaPowerTrlSimIn(:,:,outSeqCorrectLog), lowGammaPowerTrlSimOut(:,:,outSeqCorrectLog),...
    lowGammaPowerTrlSimIn(:,:,outSeqIncorrectLog), lowGammaPowerTrlSimOut(:,:,outSeqIncorrectLog), 'LowGamma-Power');
PlotMeanPowerTrlSimilarityMatrix(lowGammaPhaseTrlSimIn(:,:,inSeqCorrectLog), lowGammaPhaseTrlSimOut(:,:,inSeqCorrectLog),...
    lowGammaPhaseTrlSimIn(:,:,inSeqIncorrectLog), lowGammaPhaseTrlSimOut(:,:,inSeqIncorrectLog),...
    lowGammaPhaseTrlSimIn(:,:,outSeqCorrectLog), lowGammaPhaseTrlSimOut(:,:,outSeqCorrectLog),...
    lowGammaPhaseTrlSimIn(:,:,outSeqIncorrectLog), lowGammaPhaseTrlSimOut(:,:,outSeqIncorrectLog), 'LowGamma-Phase');

% Plot High Gamma Power Similarity Matrices
PlotMeanPowerTrlSimilarityMatrix(highGammaPowerTrlSimIn(:,:,inSeqCorrectLog), highGammaPowerTrlSimOut(:,:,inSeqCorrectLog),...
    highGammaPowerTrlSimIn(:,:,inSeqIncorrectLog), highGammaPowerTrlSimOut(:,:,inSeqIncorrectLog),...
    highGammaPowerTrlSimIn(:,:,outSeqCorrectLog), highGammaPowerTrlSimOut(:,:,outSeqCorrectLog),...
    highGammaPowerTrlSimIn(:,:,outSeqIncorrectLog), highGammaPowerTrlSimOut(:,:,outSeqIncorrectLog), 'HighGamma-Power');
PlotMeanPowerTrlSimilarityMatrix(highGammaPhaseTrlSimIn(:,:,inSeqCorrectLog), highGammaPhaseTrlSimOut(:,:,inSeqCorrectLog),...
    highGammaPhaseTrlSimIn(:,:,inSeqIncorrectLog), highGammaPhaseTrlSimOut(:,:,inSeqIncorrectLog),...
    highGammaPhaseTrlSimIn(:,:,outSeqCorrectLog), highGammaPhaseTrlSimOut(:,:,outSeqCorrectLog),...
    highGammaPhaseTrlSimIn(:,:,outSeqIncorrectLog), highGammaPhaseTrlSimOut(:,:,outSeqIncorrectLog), 'HighGamma-Phase');

% Plot Ripple Power Similarity Matrices
PlotMeanPowerTrlSimilarityMatrix(ripplePowerTrlSimIn(:,:,inSeqCorrectLog), ripplePowerTrlSimOut(:,:,inSeqCorrectLog),...
    ripplePowerTrlSimIn(:,:,inSeqIncorrectLog), ripplePowerTrlSimOut(:,:,inSeqIncorrectLog),...
    ripplePowerTrlSimIn(:,:,outSeqCorrectLog), ripplePowerTrlSimOut(:,:,outSeqCorrectLog),...
    ripplePowerTrlSimIn(:,:,outSeqIncorrectLog), ripplePowerTrlSimOut(:,:,outSeqIncorrectLog), 'Ripple-Power');
PlotMeanPowerTrlSimilarityMatrix(ripplePhaseTrlSimIn(:,:,inSeqCorrectLog), ripplePhaseTrlSimOut(:,:,inSeqCorrectLog),...
    ripplePhaseTrlSimIn(:,:,inSeqIncorrectLog), ripplePhaseTrlSimOut(:,:,inSeqIncorrectLog),...
    ripplePhaseTrlSimIn(:,:,outSeqCorrectLog), ripplePhaseTrlSimOut(:,:,outSeqCorrectLog),...
    ripplePhaseTrlSimIn(:,:,outSeqIncorrectLog), ripplePhaseTrlSimOut(:,:,outSeqIncorrectLog), 'Ripple-Phase');


fprintf('\n\n Analysis Complete      '); 
toc(ticMain)

%% Mean Similarity Plotting Function
function PlotMeanPowerTrlSimilarityMatrix(iscIN, iscOUT, isiIN, isiOUT, oscIN, oscOUT, osiIN, osiOUT, bandName)
    figure;
    subplot(2,4,1)
    imagesc(mean(iscIN,3), [-1 1]);
    hold on
    AddLines(size(iscIN,1));
    title([bandName ' InSeq Correct: Poke In']);
    subplot(2,4,2)
    imagesc(mean(iscOUT,3), [-1 1]);
    hold on
    AddLines(size(iscOUT,1));
    title('InSeq Correct: Poke Out');
    subplot(2,4,3)
    imagesc(mean(isiIN,3), [-1 1]);
    hold on
    AddLines(size(isiIN,1));
    title('InSeq Incorrect: Poke In');
    subplot(2,4,4)
    imagesc(mean(isiOUT,3), [-1 1]);
    hold on
    AddLines(size(isiOUT,1));
    title('InSeq Incorrect: Poke Out');
    subplot(2,4,5)
    imagesc(mean(oscIN,3), [-1 1]);
    hold on
    AddLines(size(oscIN,1));
    title('OutSeq Correct: Poke In');
    subplot(2,4,6)
    imagesc(mean(oscOUT,3), [-1 1]);
    hold on
    AddLines(size(oscOUT,1));
    title('OutSeq Correct: Poke Out');
    subplot(2,4,7)
    imagesc(mean(osiIN,3), [-1 1]);
    hold on
    AddLines(size(osiIN,1));
    title('OutSeq Incorrect: Poke In');
    subplot(2,4,8)
    imagesc(mean(osiOUT,3), [-1 1]);
    hold on
    AddLines(size(osiOUT,1));
    title('OutSeq Incorrect: Poke Out');
end

function AddLines(numChans)
    for lineNum = 0:4:numChans
        line([lineNum+0.5 lineNum+0.5], [0.5 numChans+0.5], 'linewidth', 2, 'color', 'white');
        line([0.5 numChans+0.5], [lineNum+0.5 lineNum+0.5], 'linewidth', 2, 'color', 'white');
    end
end