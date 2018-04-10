function FratioAnalysis_SM(unitID, grouping1, grouping2)
%% FratioAnalysis_SM.m
%
%   Comparing F-Ratio to assess information content... better description
%   needed later.
%
%%
groupType1IDs = fieldnames(grouping1);
groupType2IDs = fieldnames(grouping2);

if ~(length(groupType1IDs) == length(groupType2IDs))
    error('Groupings unequal for %s', unitID);
end

%%
% for grp = 1:length(groupType1IDs)
%     curGrp1 = grouping1.(groupType1IDs{grp});
%     curGrp2 = grouping2.(groupType2IDs{grp});