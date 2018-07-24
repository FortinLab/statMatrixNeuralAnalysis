close all
clear all
%%
origDir = cd;
[fileDir] = uigetdir(origDir);
if fileDir==0
    disp('Analysis Cancelled')
    return
else
    cd(fileDir)
end

%%
dirContents = dir(fileDir);
fileNames = {dirContents.name};

load(fileNames{cellfun(@(a)~isempty(a), strfind(fileNames, 'BehaviorMatrix'))});
load(fileNames{cellfun(@(a)~isempty(a), strfind(fileNames, 'EnsembleMatrix'))});
ensembleUniDta = ensembleMatrix(:,2:end);
samp = mode(diff(ensembleMatrix(:,1)));
%% Extract Spike Data Based on Trial Period %%
pokeInAlignedBehavMatrix = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [0 0], 'PokeIn');
% trialEnsembles = ExtractTrialData_SM(pokeInAlignedBehavMatrix, ensembleMatrix(:,2:end));

pokeOutAlignedBehavMatrix = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [0 0], 'PokeOut');
% trialEnsembles = ExtractTrialData_SM(pokeOutAlignedBehavMatrix, ensembleMatrix(:,2:end));

trialBehavMatrix = CombineTrialLogVects_SM(pokeInAlignedBehavMatrix, pokeOutAlignedBehavMatrix);
trialEnsembles = ExtractTrialData_SM(trialBehavMatrix, ensembleUniDta);

%% Create Logical Vectors & Runtime Variables
odorAlog = [trialBehavMatrix.Odor] == 1;
odorBlog = [trialBehavMatrix.Odor] == 2;
odorClog = [trialBehavMatrix.Odor] == 3;
odorDlog = [trialBehavMatrix.Odor] == 4;
odorElog = [trialBehavMatrix.Odor] == 5;

perfLog = [trialBehavMatrix.Performance] == 1;
inSeqLog = [trialBehavMatrix.ItemItemDistance] == 1;

xCorrBins = -0.2:0.01:0.2;
%% Straight xCorr Per Session
xCorr = cell(size(ensembleUniDta,2));
for u1 = 1:size(ensembleUniDta,2)
    fprintf('Unit %s Leading\n', ensembleMatrixColIDs{u1+1});
    for u2 = 1:size(ensembleUniDta,2)
        if u2>=u1
            u1SpkTms = find(ensembleUniDta(:,u1)>=1);
            u2SpkTms = find(ensembleUniDta(:,u2)>=1);
            u2RelSpkTms = nan(length(u1SpkTms),2);
            for spk1 = 2:length(u1SpkTms)-1
                postSpkDiff = (u2SpkTms(find(u2SpkTms>u1SpkTms(spk1),1,'first')) - u1SpkTms(spk1))*samp;
                if ~isempty(postSpkDiff) && (u1SpkTms(spk1+1)-u1SpkTms(spk1))*samp>postSpkDiff
                    u2RelSpkTms(spk1,1) = postSpkDiff;
                end
                preSpkDiff = (u2SpkTms(find(u2SpkTms<u1SpkTms(spk1),1,'last')) - u1SpkTms(spk1))*samp;
                if ~isempty(preSpkDiff) && (u1SpkTms(spk1-1)-u1SpkTms(spk1))*samp<preSpkDiff
                    u2RelSpkTms(spk1,2) = preSpkDiff;
                end
            end
            xCorr{u1,u2} = u2RelSpkTms(:);
        else
            continue
        end
    end
end
%% Trial Period xCorr
spkTmDiffTrl = repmat({cell(size(ensembleUniDta,2))}, [1 length(trialEnsembles)]);
for trl = 1:length(trialEnsembles)
    fprintf('Trial #%i\n', trl);
    curTrlEnsmbl = trialEnsembles{trl};
    for u1 = 1:size(curTrlEnsmbl,2)
        for u2 = 1:size(curTrlEnsmbl,2)
            if u2>=u1
                u1SpkTms = find(curTrlEnsmbl(:,u1)>=1);
                u2SpkTms = find(curTrlEnsmbl(:,u2)>=1);
                u2RelSpkTms = nan(length(u1SpkTms),2);
                if ~isempty(u1SpkTms) && ~isempty(u2SpkTms)
                    for spk1 = 1:length(u1SpkTms)
                        if spk1<=length(u1SpkTms)-1
                            postSpkDiff = (u2SpkTms(find(u2SpkTms>u1SpkTms(spk1),1,'first')) - u1SpkTms(spk1))*samp;
                            if ~isempty(postSpkDiff) && ((u1SpkTms(spk1+1)-u1SpkTms(spk1))*samp)>postSpkDiff
                                u2RelSpkTms(spk1,1) = postSpkDiff;
                            end
                        end
                        if spk1 >=2
                            preSpkDiff = (u2SpkTms(find(u2SpkTms<u1SpkTms(spk1),1,'last')) - u1SpkTms(spk1))*samp;
                            if ~isempty(preSpkDiff) && (u1SpkTms(spk1-1)-u1SpkTms(spk1))*samp<preSpkDiff
                                u2RelSpkTms(spk1,2) = preSpkDiff;
                            end
                        end
                    end
                end
                spkTmDiffTrl{trl}{u1,u2} = u2RelSpkTms;
            else
                continue
            end
        end
    end
end
%% Plot Overall and Odor-wise xCorr
for u1 = 1:size(ensembleUniDta,2)
    fprintf('Unit %s Leading\n', ensembleMatrixColIDs{u1+1});
    for u2 = 1:size(ensembleUniDta,2)
        if u1<u2
            fig = figure;
            % Odor A
%             curPairOAlag = cell2mat(cellfun(@(a)a{u1,u2}, spkTmDiffTrl(odorAlog & inSeqLog & perfLog), 'uniformoutput', 0)');
%             curPairOAlag(isnan(curPairOAlag)) = [];
            oA = subplot(5,4,1:2);
%             histogram(curPairOAlag, xCorrBins);
            aLagMean = mean(cell2mat(cellfun(@(a)histcounts(a, xCorrBins), cellfun(@(a)a{u1,u2}, spkTmDiffTrl(odorAlog & inSeqLog & perfLog), 'uniformoutput', 0), 'uniformoutput', 0)'));
            bar(xCorrBins(2:end)-(mode(diff(xCorrBins))/2),aLagMean,1, 'k');
            if sum(aLagMean) == 0
                aYmax = 0.000001;
            else
                aYmax = max(get(oA, 'ylim'));
            end
            title('Odor A');
            % Odor B
%             curPairOBlag = cell2mat(cellfun(@(b)b{u1,u2}, spkTmDiffTrl(odorBlog & inSeqLog & perfLog), 'uniformoutput', 0)');
%             curPairOBlag(isnan(curPairOBlag)) = [];
            oB = subplot(5,4,5:6);
%             histogram(curPairOBlag, xCorrBins);
            bLagMean = mean(cell2mat(cellfun(@(a)histcounts(a, xCorrBins), cellfun(@(b)b{u1,u2}, spkTmDiffTrl(odorBlog & inSeqLog & perfLog), 'uniformoutput', 0), 'uniformoutput', 0)'));
            bar(xCorrBins(2:end)-(mode(diff(xCorrBins))/2),bLagMean,1, 'k');
            if sum(bLagMean) == 0
                bYmax = 0.000001;
            else
                bYmax = max(get(oB, 'ylim'));
            end
            title('Odor B');
            % Odor C
%             curPairOClag = cell2mat(cellfun(@(c)c{u1,u2}, spkTmDiffTrl(odorClog & inSeqLog & perfLog), 'uniformoutput',0)');
%             curPairOClag(isnan(curPairOClag)) = [];
            oC = subplot(5,4,9:10);
%             histogram(curPairOClag, xCorrBins);
            cLagMean = mean(cell2mat(cellfun(@(a)histcounts(a, xCorrBins), cellfun(@(c)c{u1,u2}, spkTmDiffTrl(odorClog & inSeqLog & perfLog), 'uniformoutput',0), 'uniformoutput', 0)'));
            bar(xCorrBins(2:end)-(mode(diff(xCorrBins))/2),cLagMean,1, 'k');
            if sum(cLagMean) == 0
                cYmax = 0.000001;
            else
                cYmax = max(get(oC, 'ylim'));
            end
            title('Odor C');
            % Odor D
%             curPairODlag = cell2mat(cellfun(@(d)d{u1,u2}, spkTmDiffTrl(odorDlog & inSeqLog & perfLog), 'uniformoutput',0)');
%             curPairODlag(isnan(curPairODlag)) = [];
            oD = subplot(5,4,13:14);
%             histogram(curPairODlag, xCorrBins);
            dLagMean = mean(cell2mat(cellfun(@(a)histcounts(a, xCorrBins), cellfun(@(d)d{u1,u2}, spkTmDiffTrl(odorDlog & inSeqLog & perfLog), 'uniformoutput',0), 'uniformoutput', 0)'));
            bar(xCorrBins(2:end)-(mode(diff(xCorrBins))/2),dLagMean,1, 'k');
            if sum(dLagMean) == 0
                dYmax = 0.000001;
            else
                dYmax = max(get(oD, 'ylim'));
            end
            title('Odor D');
%             % Odor E
% %             curPairOElag = cell2mat(cellfun(@(e)e{u1,u2}, spkTmDiffTrl(odorElog & inSeqLog & perfLog), 'uniformoutput', 0)');
% %             curPairOElag(isnan(curPairOElag)) = [];
%             oE = subplot(5,4,17:18);
% %             histogram(curPairOElag, xCorrBins);
%             eLagMean = mean(cell2mat(cellfun(@(a)histcounts(a, xCorrBins), cellfun(@(e)e{u1,u2}, spkTmDiffTrl(odorElog & inSeqLog & perfLog), 'uniformoutput', 0), 'uniformoutput', 0)'));
%             bar(xCorrBins(2:end)-(mode(diff(xCorrBins))/2),eLagMean,1, 'k');
%             if sum(eLagMean) == 0
%                 eYmax = 0.000001;
%             else
%                 eYmax = max(get(oE, 'ylim'));
%             end
%             title('Odor E');
%             set([oA,oB,oC,oD,oE], 'ylim', [0 max([aYmax, bYmax, cYmax, dYmax, eYmax])], 'xlim', [min(xCorrBins) max(xCorrBins)]);
            set([oA,oB,oC,oD], 'ylim', [0 max([aYmax, bYmax, cYmax, dYmax])], 'xlim', [min(xCorrBins) max(xCorrBins)]);
            
            % Overall
            subplot(5,4,[3:4, 7:8])
            ovr = histogram(xCorr{u1,u2}, xCorrBins);            
            ovr.FaceColor = [1 1 1];
            title(sprintf('%s vs %s', ensembleMatrixColIDs{u1+1}, ensembleMatrixColIDs{u2+1}));
            
            % InSeq
            isClag = cell2mat(cellfun(@(isc)isc{u1,u2}, spkTmDiffTrl(inSeqLog & ~odorAlog & perfLog), 'uniformoutput', 0)');
            isClag(isnan(isClag)) = [];
            isc = subplot(5,4,15:16);
            iscH = histogram(isClag, xCorrBins);
            iscH.FaceColor = [1 1 1];
            set(gca, 'ylim', [0 max(get(gca, 'ylim'))]);
            title('InSeq Correct (sans A1)');
            
            % OutSeq
            osClag = cell2mat(cellfun(@(osc)osc{u1,u2}, spkTmDiffTrl(~inSeqLog & perfLog), 'uniformoutput', 0)');
            osClag(isnan(osClag)) = [];
            osc = subplot(5,4,19:20);
            oscH = histogram(osClag, xCorrBins);
            oscH.FaceColor = [1 1 1];
            set(gca, 'ylim', [0 max(get(gca, 'ylim'))]);
            title('OutSeq Correct');
            
%             linkaxes([isc, osc], 'xy');           
            saveas(fig, sprintf('%s_%s_10ms.png', ensembleMatrixColIDs{u1+1}, ensembleMatrixColIDs{u2+1}));
            close all
        end
    end
end
%% Xcorr Analyses...
%   1) Straight xcorr per session
%   2) Xcorr per item
%   3) Xcorr by phase
%   4) Xcorr by item & phase
%
%   how to assess 'significance' or at least 'relevance'... 
%   Overall distribution of spiking/relative spiking... simulate spiking?
%   Calculate liklihood of two neurons firing at their observed rates,
%   overlapping within a certain time bin with a particular frequency...
%   something like that... Chi^2 & K-S tests probably will be useful
%   here... K-S would work if I'm asking something about overlapping timing
%   of spiking events (since the dist cumsum will be over time)... Chi^2
%   may be ok for just a general dispersion test of things... maybe?
