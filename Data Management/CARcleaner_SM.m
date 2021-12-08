function CARcleaner_SM
smPath = uigetdir;
cd(smPath);
files = dir(smPath);
fileNames = {files.name};
smFileList = fileNames(cellfun(@(a)~isempty(a), regexp(fileNames, '_SM\>')))';
lfps = cell(size(smFileList))';
for f = 1:length(smFileList)
    temp = load(smFileList{f});
    lfps{f} = temp.statMatrix(:,2);
end
lfps = cell2mat(lfps);
gMn = mean(lfps,2);

%%
for f = 1:length(smFileList)
    temp = load(smFileList{f});
    statMatrix = [temp.statMatrix(:,1), temp.statMatrix(:,2)-gMn];
    statMatrixColIDs = temp.statMatrixColIDs(1:2);
    fName = [smFileList{f}(1:end-6) 'CARcleaned_SM.mat'];
    save(fName, 'statMatrix', 'statMatrixColIDs');
    fprintf('%s Saved\n', fName);
end
fprintf('All files saved\n');

