
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
load(fileNames{cellfun(@(a)~isempty(a), strfind(fileNames, 'BehaviorMatrix'))});
cd(fileDir);

%% Parameters
w = gausswin(51);
w = w/sum(w);
window = 0.125;
%%
for tet = 1:length(tetFiles)
    load(tetFiles{tet});    
    samp = 1/mode(diff(statMatrix(:,1)));
    windowNdx = round(window*samp);
    [bs, ~] = regexp(statMatrixColIDs, 'LFP_Raw$');
    lfpRaw = statMatrix(:,cellfun(@(a)~isempty(a), bs));
    [thetaHilb, theta] = PhaseFreqDetectAbbr(lfpRaw, statMatrix(:,1), 4, 8);
    [betaHilb, beta] = PhaseFreqDetectAbbr(lfpRaw, statMatrix(:,1), 16, 32);
%     betaRMS = zscore(conv(sqrt(conv(beta.^2, ones(windowNdx,1)/windowNdx, 'same')), w, 'same'));
    betaRMS = conv(sqrt(conv(beta.^2, ones(windowNdx,1)/windowNdx, 'same')), w, 'same');
    [betaRMShilb, ~] = PhaseFreqDetectAbbr(betaRMS, statMatrix(:,1));    
    
    pokeInTimesMtx = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [0 0], 'PokeIn');
    pokeInTimes = ExtractTrialData_SM(pokeInTimesMtx, behavMatrix(:,1));
%     pokeOutTimes = cellfun(@(a){a+0.5}, pokeInTimes);
    pokeOutTimesMtx = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [0 0], 'PokeOut');
    pokeOutTimes = ExtractTrialData_SM(pokeOutTimesMtx, behavMatrix(:,1));
%     pokeInTimes = cellfun(@(a){a-0.5}, pokeOutTimes);
    
    odorVals = [pokeInTimesMtx.Odor];
    posVals = [pokeInTimesMtx.Position];
    inSeqLog = [pokeInTimesMtx.ItemItemDistance]==1;
    performanceLog = [pokeInTimesMtx.Performance]==1;
    
    tbPwrCorr = nan(1,length(pokeInTimes));
    for trl = 1:length(pokeInTimes)
        tempPokeInNdx = find(statMatrix(:,1)==pokeInTimes{trl});
        tempPokeOutNdx = find(statMatrix(:,1)==pokeOutTimes{trl});
%         [tbPwrCorr(trl),~] = circ_corrcl(thetaHilb(tempPokeInNdx:tempPokeOutNdx), betaRMS(tempPokeInNdx:tempPokeOutNdx));
        [tbPwrCorr(trl),~] = circ_corrcc(thetaHilb(tempPokeInNdx:tempPokeOutNdx), betaRMShilb(tempPokeInNdx:tempPokeOutNdx));        
    end
    figure;
    subplot(3,1,1);
    histogram(tbPwrCorr(performanceLog), 0:0.1:1);
    title('All Trials')
    subplot(3,1,2);
    histogram(tbPwrCorr(performanceLog & inSeqLog & ~(odorVals==1)), 0:0.1:1);
    title('InSeq Correct Trials (sans odorA)');
    subplot(3,1,3);
    histogram(tbPwrCorr(performanceLog & ~inSeqLog), 0:0.1:1);
    title('OutSeq Correct Trials');
    
    figTitle = annotation('textbox', 'position', [0.025 0.935 0.7 0.05], 'String', tetFiles{tet},...
        'linestyle', 'none', 'horizontalalignment', 'left', 'interpreter', 'none');
    drawnow
end
