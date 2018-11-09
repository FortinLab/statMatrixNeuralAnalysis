%% rms calculation of an input signal determined for a given frequency. using the same calculation as in Simulink block.
function [U_rms,time_rms] = mvrms_mfile_simulink(t,y,f1,figures)
% Use this function to calculate the rms value of the input signal ('y')
% determined for a given period of time ('f1' in Hz).
% The used inputs of this function:
% 1- t: the time signal of the selected signal (in which the RMS to be
% calculated).
% 2- y: the value signal of the selected signal (as vector).
% 3- f1: the frequency in which the RMS value to be determined. is input
% will be changed to time in S insided this code.
% 4- figures: this a (0-1) input if 1 will show the results of the diffrent
% calculation steps of this function.
% The outputs of this function are:
% 1- U_rms = which is the rms values of the input Signal (as in Simulink
% the first rms valeus can be ignored).
% 2- time_rms: is the time signal of the rms values signal which as in the
% Simulink is the same as the input time signal.
[r,~] = size(y);
if r == 1
    y = y';
    t = t';
end
delta = mean(diff(t));%t(2) - t(1);
L1 = round((1/f1)/ delta);
Zeros = 0.*(0:(L1-1))';
Signal = y.*conj(y);
if figures
    figure;
    subplot(8,1,1);
    plot(t,Signal);
end
Signal = Signal * f1;
if figures
    subplot(8,1,2);
    plot(t,Signal);
end
Signal = cumtrapz(t,Signal);
if figures
    subplot(8,1,3);
    plot(t,Signal);
end
Signal1 = [Zeros;Signal(1:(length(Signal)-L1))];
if figures
    subplot(8,1,4);
    plot(t,Signal);
end
Signal = Signal - Signal1;
if figures
    subplot(8,1,5);
    plot(t,Signal);
end
Signal = abs(Signal);
if figures
    subplot(8,1,6);
    plot(t,Signal);
end
U_rms = sqrt(Signal);
if figures
    subplot(8,1,7);
    plot(t,Signal);
end
time_rms = t;
if figures
    subplot(8,1,8);
    plot(t,Signal);
end
end