function PlotPowerPhaseCorrelations_OE(iscINpwr, iscINphs,...
    iscOUTpwr, iscOUTphs,...
    isiINpwr, isiINphs,...
    isiOUTpwr, isiOUTphs,...
    oscINpwr, oscINphs,...
    oscOUTpwr, oscOUTphs,...
    osiINpwr, osiINphs,...
    osiOUTpwr, osiOUTphs,...
    bandName)
%%
figure;
%%  Comment in this section if you want to use the absolute value of rho rather than raw
% iscINpwr = abs(iscINpwr);
% iscINphs = abs(iscINphs);
% iscOUTpwr = abs(iscOUTpwr);
% iscOUTphs = abs(iscOUTphs);
% isiINpwr = abs(isiINpwr);
% isiINphs = abs(isiINphs);
% isiOUTpwr = abs(isiOUTpwr);
% isiOUTphs = abs(isiOUTphs);
% oscINpwr = abs(oscINpwr);
% oscINphs = abs(oscINphs);
% oscOUTpwr = abs(oscOUTpwr);
% oscOUTphs = abs(oscOUTphs);
% osiINpwr = abs(osiINpwr);
% osiINphs = abs(osiINphs);
% osiOUTpwr = abs(osiOUTpwr);
% osiOUTphs = abs(osiOUTphs);
%% Comment this section in if you want to run the analysis on the mean rho rather than by trial
iscINpwr = mean(iscINpwr,3);
iscINphs = mean(iscINphs,3);
iscOUTpwr = mean(iscOUTpwr,3);
iscOUTphs = mean(iscOUTphs,3);
isiINpwr = mean(isiINpwr,3);
isiINphs = mean(isiINphs,3);
isiOUTpwr = mean(isiOUTpwr,3);
isiOUTphs = mean(isiOUTphs,3);
oscINpwr = mean(oscINpwr,3);
oscINphs = mean(oscINphs,3);
oscOUTpwr = mean(oscOUTpwr,3);
oscOUTphs = mean(oscOUTphs,3);
osiINpwr = mean(osiINpwr,3);
osiINphs = mean(osiINphs,3);
osiOUTpwr = mean(osiOUTpwr,3);
osiOUTphs = mean(osiOUTphs,3);
%%
% Create Logical Matrices to remove diagonal from the plots
nonSameID = ~eye(size(iscINpwr,1));
% Remove the indices that are from the same region
for chan = 1:4:size(iscINpwr,1)
    nonSameID(chan:chan+3, chan:chan+3) = false;
end
iscLog = repmat(nonSameID, [1,1,size(iscINpwr,3)]);
oscLog = repmat(nonSameID, [1,1,size(oscINpwr,3)]);
isiLog = repmat(nonSameID, [1,1,size(isiINpwr,3)]);
osiLog = repmat(nonSameID, [1,1,size(osiINpwr,3)]);

%% Plot everything
iscIN = subplot(2,4,1);
corrScatPlot(iscINpwr(iscLog), iscINphs(iscLog), {[bandName ' Power Coupling']}, {[bandName ' Phase Coherence']},[]);
title([bandName ' InSeq Correct: Poke In']);
iscOUT = subplot(2,4,2);
corrScatPlot(iscOUTpwr(iscLog), iscOUTphs(iscLog), {[bandName ' Power Coupling']}, {[bandName ' Phase Coherence']},[]);
title('InSeq Correct: Poke Out');
isiIN = subplot(2,4,3);
corrScatPlot(isiINpwr(isiLog), isiINphs(isiLog), {[bandName ' Power Coupling']}, {[bandName ' Phase Coherence']},[]);
title('InSeq Incorrect: Poke In');
isiOUT = subplot(2,4,4);
corrScatPlot(isiOUTpwr(isiLog), isiOUTphs(isiLog), {[bandName ' Power Coupling']}, {[bandName ' Phase Coherence']},[]);
title('InSeq Incorrect: Poke Out');

oscIN = subplot(2,4,5);
corrScatPlot(oscINpwr(oscLog), oscINphs(oscLog), {[bandName ' Power Coupling']}, {[bandName ' Phase Coherence']},[]);
title('OutSeq Correct: Poke In');
oscOUT = subplot(2,4,6);
corrScatPlot(oscOUTpwr(oscLog), oscOUTphs(oscLog), {[bandName ' Power Coupling']}, {[bandName ' Phase Coherence']},[]);
title('OutSeq Correct: Poke Out');
osiIN = subplot(2,4,7);
corrScatPlot(osiINpwr(osiLog), osiINphs(osiLog), {[bandName ' Power Coupling']}, {[bandName ' Phase Coherence']},[]);
title('OutSeq Incorrect: Poke In');
osiOUT = subplot(2,4,8);
corrScatPlot(osiOUTpwr(osiLog), osiOUTphs(osiLog), {[bandName ' Power Coupling']}, {[bandName ' Phase Coherence']},[]);
title('OutSeq Incorrect: Poke OUT');

linkaxes([iscIN, iscOUT, isiIN, isiOUT, oscIN, oscOUT, osiIN, osiOUT], 'xy');

end


