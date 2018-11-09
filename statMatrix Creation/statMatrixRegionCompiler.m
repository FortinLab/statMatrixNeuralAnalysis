function statMatrixRegionCompiler(fileDir,targetFolder)
% fileDir = location of individual channel statMatrix files
% targetFolder = location of resulting compiled statMatrices
%% Detects the number of statMatrix files in the fileDir
origDir = cd; % saves current file directory
cd(fileDir) % changes to target directory
FileName = 'statMatrix_';
numfiles = length(dir([FileName,'*']));
cd(origDir) % reverts to original directory

for channelNum = 1:numfiles
    %% Load current statMatrix file onto workspace
    currFileName = sprintf('statMatrix_%d',channelNum);
    statMatrixData = load(currFileName);
    statMatrix = statMatrixData.statMatrix;
    statMatrixCols = statMatrixData.statMatrixColIDs;
    
    %% Determine how many LFP and Behav columns are there saves some primary files
    % Since the columns containing behavioral events are binary, we can use
    % that to differentiate between them and the LFP columns
    if channelNum == 1
        zeroCountCol = sum(~statMatrix);
        numBehavCol = length(find(zeroCountCol>1));
        numLFPCol = size(statMatrix,2)-numBehavCol;
        
        %Saves column labels
        Behav_Matrix_ColIDs = statMatrixCols(:,numLFPCol+1:length(zeroCountCol));
        BehavLabelfilepath = strcat(targetFolder,'\','Behav_Matrix_ColIDs','.mat');
        save(BehavLabelfilepath,'Behav_Matrix_ColIDs')
        
        LFP_Matrix_ColIDs = statMatrixCols(:,1:numLFPCol);
        LFPLabelfilepath = strcat(targetFolder,'\','LFP_Matrix_ColIDs','.mat');
        save(LFPLabelfilepath,'LFP_Matrix_ColIDs')
        
        
        % Saves the behavioral columns file in this conditional, so that it
        % only happens once
        BehavEvents = statMatrix(:,numLFPCol+1:length(zeroCountCol));
        Behavfilepath = strcat(targetFolder,'\','BehavEvents','.mat');
        save(Behavfilepath,'BehavEvents')
    end
    
    %% Stores LFP arrays onto a sections of 16
    % Did this to avoid an array memory error for my computer when it
    % maximizes in RAM. Looks pretty disgusting, as I'm not sure how to
    % deal with this other than storing all statMatrices into a single
    % array
    
    statMatrix = statMatrix(:,1:numLFPCol);
    if channelNum <= 16
        statMatrixAll_A(:,:,channelNum) = statMatrix;
    elseif channelNum >16 && channelNum <=32
        statMatrixAll_B(:,:,channelNum-16) = statMatrix;
    elseif channelNum >32 && channelNum <=48
        statMatrixAll_C(:,:,channelNum-32) = statMatrix;
    elseif channelNum >48 && channelNum <=64
        statMatrixAll_D(:,:,channelNum-48) = statMatrix;
    end
    
    if mod(channelNum,4) == 0
        if channelNum <= 16
            statMatrix_Region = statMatrixAll_A(:,:,channelNum-3:channelNum);
        elseif channelNum >16 && channelNum <=32
            statMatrix_Region = statMatrixAll_B(:,:,channelNum-19:channelNum-16);
        elseif channelNum >32 && channelNum <=48
            statMatrix_Region = statMatrixAll_C(:,:,channelNum-35:channelNum-32);
        elseif channelNum >48 && channelNum <=64
            statMatrix_Region = statMatrixAll_D(:,:,channelNum-51:channelNum-48);
        end
        filename = sprintf('statMatrixRegion_%d',channelNum/4);
        fullfilepath = strcat(targetFolder,'\',filename,'.mat');
        save(fullfilepath,'statMatrix_Region')
    end
    
end

end
