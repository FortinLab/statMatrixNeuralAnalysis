

order = 3;
fcut  = 60; 
NyquistFreq = 1000/2;
[b,l] = butter(order,fcut/(NyquistFreq),'low'); 
LFP_Bandpass = filtfilt(b,l,LFP_Raw);