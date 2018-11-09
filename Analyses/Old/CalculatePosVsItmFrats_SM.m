function eveFvals = CalculatePosVsItmFrats_SM(uniSpkData, behavMatrixTrls, eventWindows, pehBinSize, fMax)
%%
%
%
%%
odorIDs = [behavMatrixTrls.Odor];
posIDs = [behavMatrixTrls.Position];

eveFvals = cell(size(uniSpkData,2), 3);
for u = 1:size(uniSpkData,2)
    tempOdorIDs = odorIDs;
    tempPosIDs = posIDs;
    curUniData = ExtractTrialData_SM(behavMatrixTrls, uniSpkData(:,u));
    % Remove Empty Trials
    noSpkLog = cellfun(@(a)isempty(a), curUniData);
    curUniData(noSpkLog) = [];
    tempOdorIDs(noSpkLog) = [];
    tempPosIDs(noSpkLog) = [];
    curUniPEH = cell(length(curUniData),1);
    for trl = 1:length(curUniPEH)
        [curUniPEH{trl}, ~] = RebinPEH_SM(curUniData{trl}, eventWindows, pehBinSize);
    end
    curUniUnpacked = cell2mat(curUniPEH);
    odorF = nan(1,size(curUniUnpacked,2));
    positionF = nan(1,size(curUniUnpacked,2));
    for ndx = 1:size(curUniUnpacked,2)
        [~,odrTable,~] = anova1(curUniUnpacked(:,ndx), tempOdorIDs', 'off');
        curOdorF = odrTable{2,5};
        if ~isempty(curOdorF) && curOdorF>0 && curOdorF<fMax
            odorF(ndx) = curOdorF;
        elseif curOdorF>fMax
            odorF(ndx) = fMax;
        else
            odorF(ndx) = nan;
        end
        [~,posTable,~] = anova1(curUniUnpacked(:,ndx), tempPosIDs', 'off');
        curPosF = posTable{2,5};
        if ~isempty(curPosF) && curPosF>0 && curPosF<fMax
            positionF(ndx) = curPosF;
        elseif curPosF>fMax
            positionF(ndx) = fMax;
        else
            positionF(ndx) = nan;
        end
    end
    eveFvals(u,:) = [{positionF}, {odorF}, {positionF-odorF}];
end