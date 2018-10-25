function [rawFvals, chanceFvals, zNormFvals] = UnitFvalCalcPERM_SM(unitSMmatrix, idVect, dataBinSize, numPerms)
%% UnitFvalCalcPERM
%   Uses a sliding window analysis to evaluate how well the activity for
%   each unit is accounted for by the identity of the trial it was recorded
%   during. It first calculates the observed F-ratio values at each time
%   point and then 
%
%   Inputs:
%       - unitSMmatrix: unit data compiled from the statMatrix
%           organization. Sized MxNxO where M=timestamp, N=unit# and
%           O=trial#. Typically extracted using EpochExtraction_SM.
%       - idVect: vector identifying the trial IDs. Should be the same
%           length as the third dimension of unitSMmatrix.
%       - dataBinSize: length of the time period used to bin the data.
%           Future iterations may include variable inputs to augment this
%           selection, e.g. LFP phase or something like that. Default =
%           200.
%       - numPerms: number of permutations used to calculate the chance
%           distribution. Default = 100.
%
%   Outputs:
%       - rawFvals: Observed F-ratios over time for each unit organized in
%           an MxN matrix where M=timestamp and N=unitID#.
%       - chanceFvals: Chance distribution of F-ratios over time for each
%           unit organized in an MxNxO matrix where M=timestamp, N=unitID#
%           and O=permutation#.
%       - zNormFvals: Z-score of the observed F-ratios overtime relative to
%           the chance distribution.
%
%   Notes:
%   The time periods under investigation are determined by the number of
%   rows in the unitSMmatrix (number of timestamps) and the dataBinSize
%   input. Specifically, the number of rows in unitSMmatrix should
%   correspond to the time period of interest +/- (dataBinSize/2). This is
%   to ensure accurate calculations for the time periods of interest.
%       For example: 
%           To examine a period -800ms : 500ms relative to poke initiation
%           using a 200ms window, unitSMmatrix should contain data from
%           -900ms : 600ms relative to poke initiation (100ms before and
%           after the desired analysis window).
%% Chance up that random stream
randperm(floor(sum(clock)));
%% Calculate Runtime variables
timeTrimmingMask = false(size(unitSMmatrix,1),1);

if mod(dataBinSize,2)==0
    halfBin = dataBinSize/2;
    timeTrimmingMask(size(unitSMmatrix,1)-(halfBin-1):end) = true;
else
    halfBin = floor(dataBinSize/2);
    timeTrimmingMask(size(unitSMmatrix,1)-(halfBin):end) = true;
end
timeTrimmingMask(1:halfBin) = true;
%% Set up outputs
rawFvals = nan(size(unitSMmatrix,1)-(dataBinSize), size(unitSMmatrix,2));
chanceFvals = nan(size(unitSMmatrix,1)-(dataBinSize), size(unitSMmatrix,2), numPerms);
zNormFvals = nan(size(unitSMmatrix,1)-(dataBinSize), size(unitSMmatrix,2));

%% Run things
for uni = 1:size(unitSMmatrix,2)
%     tic
    % First bin the spike counts to speed up calculations below
    curUniEnsmblBinned = nan(size(unitSMmatrix,1), size(unitSMmatrix,3));
    for trl = 1:size(unitSMmatrix,3)                                                                        % Data is binned on a per trial basis
        if mod(dataBinSize,2)==0
            curUniEnsmblBinned(:,trl) = conv(unitSMmatrix(:,uni,trl)', ones(1,dataBinSize+1), 'same');        % I like to convolve with an odd number
        else
            curUniEnsmblBinned(:,trl) = conv(unitSMmatrix(:,uni,trl)', ones(1,dataBinSize), 'same');        
        end
    end
    curUniEnsmblBinned(timeTrimmingMask,:) = [];
    parfor t = 1:size(curUniEnsmblBinned,1)
        [~,table,~] = anova1(curUniEnsmblBinned(t,:)', idVect, 'off');
        if ~isempty(table{2,5})
            rawFvals(t,uni) = table{2,5};
        else
            rawFvals(t,uni) = 1;
        end
    end
    rawFvals(isnan(rawFvals(:,uni)),uni) = 1;                                                               % Replace any NaN values with 1 (chance F-ratio value)
    for r = 1:numPerms
        idVectShflVect = randperm(length(idVect));
        idVectShfl = idVect(idVectShflVect);
        parfor t = 1:size(curUniEnsmblBinned,1)
            [~,tableRND,~] = anova1(curUniEnsmblBinned(t,:)', idVectShfl, 'off');
            if ~isempty(tableRND{2,5})
                chanceFvals(t,uni,r) = tableRND{2,5};
            else
                chanceFvals(t,uni,r) = 1;
            end
        end
        chanceFvals(isnan(chanceFvals(:,uni,r)),uni,r) = 1;                                                 % Replace any NaN values with 1 (chance F-ratio value)
    end
    for t = 1:size(curUniEnsmblBinned,1)
        zUnifFvect = zscore([rawFvals(t,uni), reshape(chanceFvals(t,uni,:), [1 size(chanceFvals,3)])]);
        zNormFvals(t,uni) = zUnifFvect(1);
    end
%     toc
end

    
    
    
    
    
    
    