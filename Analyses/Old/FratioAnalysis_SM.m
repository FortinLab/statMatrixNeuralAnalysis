function [statDiff] = FratioAnalysis_SM(unitID, grouping1, grouping2)
%% FratioAnalysis_SM.m
%
%   Comparing F-Ratio to assess information content... better description
%   needed later.
%
%%

% Check first level of the grouping structures
group1EventIDs = fieldnames(grouping1);
group2EventIDs = fieldnames(grouping2);
eventIDs = unique([group1EventIDs,group2EventIDs]);
if ~(length(group1EventIDs) == length(group2EventIDs))
    error('Events unequal for %s', unitID);
elseif sum(strcmp(eventIDs,intersect(group1EventIDs, group2EventIDs)))<length(group1EventIDs)
    error('Non-identical enevts for %s', unitID);
end

% Check second level of the grouping structures
group1EventGroupingIDs = fieldnames(grouping1.(eventIDs{1}));
group2EventGroupingIDs = fieldnames(grouping2.(eventIDs{1}));
eventGroupingIDs = unique([group1EventGroupingIDs, group2EventGroupingIDs]);
if ~(length(group1EventGroupingIDs) == length(group2EventGroupingIDs))
    error('Groupings unequal for %s', unitID);
elseif sum(strcmp(eventGroupingIDs,intersect(group1EventGroupingIDs, group2EventGroupingIDs)))<length(group1EventGroupingIDs)
    error('Non-identical groupings for %s', unitID);
end

if isstruct(grouping1.(eventIDs{1}).(eventGroupingIDs{1}))
    comparison = 'F-test';
    group1TrialTypeIDs = fieldnames(grouping1.(eventIDs{1}).(eventGroupingIDs{1}));
    group2TrialTypeIDs = fieldnames(grouping2.(eventIDs{1}).(eventGroupingIDs{1}));
    
    if ~(length(group1TrialTypeIDs)==length(group2TrialTypeIDs))
        error('Trial types unequal for %s', unitID);
    end
else
    comparison = 'T-test';
end

%% Run though each event and grouping and compute statistics
for eve = 1:length(eventIDs)
    curEvent = eventIDs{eve};
    for grp = 1:length(eventGroupingIDs)
        curEveGrp = eventGroupingIDs{grp};
        switch comparison
            case 'F-test'
                grp1Data = struct2cell(grouping1.(curEvent).(curEveGrp));
                grp1TrlLogs = arrayfun(@(a,b)ones(length(a{1}),1)*b, grp1Data, (1:length(grp1Data))', 'uniformoutput', 0);
                grp1DataUnpacked = cell2mat(vertcat(grp1Data{:}));
                grp1TrlLogsUnpacked = cell2mat(grp1TrlLogs);
                grp1F = nan(1,size(grp1DataUnpacked,2));
                for t1 = 1:size(grp1DataUnpacked,2)
                    [~,grp1Table,~] = anova1(grp1DataUnpacked(:,t1), grp1TrlLogsUnpacked, 'off');
                    curF = grp1Table{2,5};
                    if ~isempty(curF)
                        grp1F(t1) = curF;
                    else
                        grp1F(t1) = 0;
                    end
                end
                grp2Data = struct2cell(grouping2.(curEvent).(curEveGrp));
                grp2TrlLogs = arrayfun(@(a,b)ones(length(a{1}),1)*b, grp2Data, (1:length(grp2Data))', 'uniformoutput', 0);
                grp2DataUnpacked = cell2mat(vertcat(grp2Data{:}));
                grp2TrlLogsUnpacked = cell2mat(grp2TrlLogs);
                grp2F = nan(1,size(grp2DataUnpacked,2));
                for t2 = 1:size(grp2DataUnpacked,2)
                    [~,grp2Table,~] = anova1(grp2DataUnpacked(:,t2), grp2TrlLogsUnpacked, 'off');
                    curF = grp2Table{2,5};
                    if ~isempty(curF)
                        grp2F(t2) = curF;
                    else
                        grp2F(t2) = 0;
                    end
                end
                statDiff.(curEvent).(curEveGrp) = grp2F-grp1F;                
            case 'T-test'
        end
    end
end
%%