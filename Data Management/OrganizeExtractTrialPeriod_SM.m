function behavTrialMtx = OrganizeExtractTrialPeriod_SM(behavMatrix, behavMatrixColIDs, statMatrix, preStartDur, postEndDur)
% Combines the Organize and Extract functions to pull out trial data from
% the statMatrix input within the window specified by preStartDur and
% postEndDur reflecting the amount of time prior to trial initiation and
% after trial end you want extracted.

sampleRate = 1/mode(diff(behavMatrix(:,1)));
trlWindow = [round(preStartDur*sampleRate) round(postEndDur*sampleRate)];

behavMatrixTrialStruct = OrganizeTrialData_SM(behavMatrix, behavMatrixColIDs, [0 0], 'PokeIn');
behavTrialMtx = rmfield(behavMatrixTrialStruct, 'TrialLogVect');
for t = 1:length(behavMatrixTrialStruct)
    behavTrialMtx(t).StatMatrixData = statMatrix(behavMatrixTrialStruct(t).PokeInIndex - trlWindow(1) : behavMatrixTrialStruct(t).PokeOutIndex + trlWindow(2),:);
end