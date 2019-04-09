function PlotPowerPhaseCorrelationsMatrix_OE(iscINpwr, iscINphs,...
    iscOUTpwr, iscOUTphs,...
    isiINpwr, isiINphs,...
    isiOUTpwr, isiOUTphs,...
    oscINpwr, oscINphs,...
    oscOUTpwr, oscOUTphs,...
    osiINpwr, osiINphs,...
    osiOUTpwr, osiOUTphs,...
    bandName, cLims)
%% Calculate Correllation Matrices
iscIN = nan(size(iscINpwr,1));
iscOUT = nan(size(iscOUTpwr,1));
isiIN = nan(size(isiINpwr,1));
isiOUT = nan(size(isiOUTpwr,1));
oscIN = nan(size(oscINpwr,1));
oscOUT = nan(size(oscOUTpwr,1));
osiIN = nan(size(osiINpwr,1));
osiOUT = nan(size(osiOUTpwr,1));

for ndx1 = 1:size(iscINpwr,1)
    for ndx2 = 1:size(iscINpwr,1)
        iscIN(ndx1,ndx2) = corr(reshape(iscINpwr(ndx1,ndx2,:), [size(iscINpwr,3),1]), reshape(iscINphs(ndx1,ndx2,:), [size(iscINphs,3),1]));
        iscOUT(ndx1,ndx2) = corr(reshape(iscOUTpwr(ndx1,ndx2,:), [size(iscOUTpwr,3),1]), reshape(iscOUTphs(ndx1,ndx2,:), [size(iscOUTphs,3),1]));
        isiIN(ndx1,ndx2) = corr(reshape(isiINpwr(ndx1,ndx2,:), [size(isiINpwr,3),1]), reshape(isiINphs(ndx1,ndx2,:), [size(isiINphs,3),1]));
        isiOUT(ndx1,ndx2) = corr(reshape(isiOUTpwr(ndx1,ndx2,:), [size(isiOUTpwr,3),1]), reshape(isiOUTphs(ndx1,ndx2,:), [size(isiOUTphs,3),1]));
        oscIN(ndx1,ndx2) = corr(reshape(oscINpwr(ndx1,ndx2,:), [size(oscINpwr,3),1]), reshape(oscINphs(ndx1,ndx2,:), [size(oscINphs,3),1]));
        oscOUT(ndx1,ndx2) = corr(reshape(oscOUTpwr(ndx1,ndx2,:), [size(oscOUTpwr,3),1]), reshape(oscOUTphs(ndx1,ndx2,:), [size(oscOUTphs,3),1]));
        osiIN(ndx1,ndx2) = corr(reshape(osiINpwr(ndx1,ndx2,:), [size(osiINpwr,3),1]), reshape(osiINphs(ndx1,ndx2,:), [size(osiINphs,3),1])); 
        osiOUT(ndx1,ndx2) = corr(reshape(osiOUTpwr(ndx1,ndx2,:), [size(osiOUTpwr,3),1]), reshape(osiOUTphs(ndx1,ndx2,:), [size(osiOUTphs,3),1]));
    end
end

%% Plot everything
figure;
subplot(2,4,1)
imagesc(iscIN, cLims);
hold on
AddLines(size(iscIN,1));
title([bandName ' InSeq Correct: Poke In']);
subplot(2,4,2)
imagesc(iscOUT, cLims);
hold on
AddLines(size(iscOUT,1));
title('InSeq Correct: Poke Out');
subplot(2,4,3)
imagesc(isiIN, cLims);
hold on
AddLines(size(isiIN,1));
title('InSeq Incorrect: Poke In');
subplot(2,4,4)
imagesc(isiOUT, cLims);
hold on
AddLines(size(isiOUT,1));
title('InSeq Incorrect: Poke Out');
subplot(2,4,5)
imagesc(oscIN, cLims);
hold on
AddLines(size(oscIN,1));
title('OutSeq Correct: Poke In');
subplot(2,4,6)
imagesc(oscOUT, cLims);
hold on
AddLines(size(oscOUT,1));
title('OutSeq Correct: Poke Out');
subplot(2,4,7)
imagesc(osiIN, cLims);
hold on
AddLines(size(osiIN,1));
title('OutSeq Incorrect: Poke In');
subplot(2,4,8)
imagesc(osiOUT, cLims);
hold on
AddLines(size(osiOUT,1));
title('OutSeq Incorrect: Poke Out');


    function AddLines(numChans)
        for lineNum = 0:4:numChans
            line([lineNum+0.5 lineNum+0.5], [0.5 numChans+0.5], 'linewidth', 1, 'color', 'white');
            line([0.5 numChans+0.5], [lineNum+0.5 lineNum+0.5], 'linewidth', 1, 'color', 'white');
        end
    end
end