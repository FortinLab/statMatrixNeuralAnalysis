function PlotMeanPowerTrlSimilarityMatrix_OE(iscIN, iscOUT, isiIN, isiOUT, oscIN, oscOUT, osiIN, osiOUT, bandName, cLims)
figure;
subplot(2,4,1)
imagesc(median(iscIN,3), cLims);
hold on
AddLines(size(iscIN,1));
title([bandName ' InSeq Correct: Poke In']);
subplot(2,4,2)
imagesc(median(iscOUT,3), cLims);
hold on
AddLines(size(iscOUT,1));
title('InSeq Correct: Poke Out');
subplot(2,4,3)
imagesc(median(isiIN,3), cLims);
hold on
AddLines(size(isiIN,1));
title('InSeq Incorrect: Poke In');
subplot(2,4,4)
imagesc(median(isiOUT,3), cLims);
hold on
AddLines(size(isiOUT,1));
title('InSeq Incorrect: Poke Out');
subplot(2,4,5)
imagesc(median(oscIN,3), cLims);
hold on
AddLines(size(oscIN,1));
title('OutSeq Correct: Poke In');
subplot(2,4,6)
imagesc(median(oscOUT,3), cLims);
hold on
AddLines(size(oscOUT,1));
title('OutSeq Correct: Poke Out');
subplot(2,4,7)
imagesc(median(osiIN,3), cLims);
hold on
AddLines(size(osiIN,1));
title('OutSeq Incorrect: Poke In');
subplot(2,4,8)
imagesc(median(osiOUT,3), cLims);
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