% Input Variables

recStart = 56576;
processorNum = 100;
channelNum = 1;
fileDir = 'C:\Users\Gabo\Documents\MATLAB\NCP04_2017-11-03_09-20-28_task217';
ssnDir = 'C:\Users\Gabo\Documents\MATLAB\Good Sessions for Code Dev\NCP04_SESSION217_03-Nov-2017';

% Load data

[~,statMatrix, ~] = statMatrixOE(recStart, processorNum ,channelNum, fileDir, ssnDir, 'processor', 'ADC', 'CH','session');
LFP_raw = statMatrix(:,2);

% Calc rms

LFP_rms = rms(LFP_raw);
plot(LFP_rms)




