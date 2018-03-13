function trialOrgData = ExtractTrialData_SM(behavMatrixTrialStruct, inputMatrix)
%% ExtractTrialData_SM
%   Extracts trial-wise data from inputMatrix using the trial logical
%   vectors contained within the behavMatrixTrialStruct structure, created
%   using OrganizeTrialData_SM.m
%
%% 
trialOrgData = cell(1,length(behavMatrixTrialStruct));
for trl = 1:length(behavMatrixTrialStruct)
    trialOrgData{trl} = inputMatrix(behavMatrixTrialStruct(trl).TrialLogVect,:);
end