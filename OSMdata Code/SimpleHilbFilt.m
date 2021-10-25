function [zeroTroughPhase, rawPow, zPow] = SimpleHilbFilt(data, samp, window)
    Wn_FRange = [window(1)/(samp/2) window(2)/(samp/2)]; % normalized by the nyquist frequency
    [bFRange, aFRange] = butter(3,Wn_FRange);
    filtered = filtfilt(bFRange,aFRange,data);

    zeroTroughPhase = atan2(imag(hilbert(filtered*-1)), real(hilbert(filtered*-1)));
    rawPow = abs(hilbert(filtered));
    zPow = zscore(rawPow);
end