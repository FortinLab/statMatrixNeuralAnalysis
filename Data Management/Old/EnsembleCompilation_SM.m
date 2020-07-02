function EnsembleCompilation_SM
%% EnsembleCompilation_SM
%   Code to extract all the unit activity from a set of files in the
%   statMatrix format FROM THE SAME SESSION and compile them into a common
%   variable
%
% 02/14/2018    - Created by GE
%
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
tetFileLog = cellfun(@(a)~isempty(a), regexp(fileNames, '_T([0-9]*).mat'));
tetFiles = fileNames(tetFileLog)';
load(tetFiles{1}, 'statMatrixColIDs');

%% Compile everything!
ensembleMatrix = cell(1,length(tetFiles));
ensembleMatrixColIDs = cell(1,length(tetFiles));
ensembleUnitSummaries = cell(1,length(tetFiles));
for tet = 1:length(tetFiles)
    load(tetFiles{tet});
    ensembleUnitSummaries{tet} = unitSummary;
    if tet == 1
        tsVect = statMatrix(:,1); %#ok<NODEF>
    end
    statMatrixColIDs = statMatrixColIDs(cellfun(@(a)~isempty(a), statMatrixColIDs));
    uniColLog = cellfun(@(a)~isempty(a), regexp(statMatrixColIDs, '-U([0-9]*)'));
    ensembleMatrix{tet} = statMatrix(:,uniColLog);
    ensembleMatrixColIDs{tet} = statMatrixColIDs(uniColLog);
    fprintf('%s Complete\n', tetFiles{tet});
end
%%
ensembleMatrix = [tsVect, cell2mat(ensembleMatrix)]; %#ok<NASGU>
ensembleMatrixColIDs = [{'TimeBin'}, ensembleMatrixColIDs{:}]; %#ok<NASGU>
ensembleUnitSummaries = cell2mat(ensembleUnitSummaries); %#ok<NASGU>

save('EnsembleMatrix.mat', 'ensembleMatrix', 'ensembleMatrixColIDs', 'ensembleUnitSummaries', '-v7.3');

    
