InSeqLog_Session = [ssnData.TranspositionDistance];
InSeqLog = zeros(1,size(CH_Activation_TS.OnPeakTS_1,2));

for k = 1:size(InSeqLog_Session,2)
    if InSeqLog_Session(k) == 0
       InSeqLog(AllOdor_Index(k)) = 1;
    else
       InSeqLog(AllOdor_Index(k)) = 0;
    end
end

InSeqLog = InSeqLog';
