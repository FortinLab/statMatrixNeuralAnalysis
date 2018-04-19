% Question 9.2
disp('Locate eog_trimmed.mat in your device')
uiopen('.mat') % Look for eog_trimmed.mat

% 9.2.a

disp('Question 9.2.a')
disp('The first figure refers to the signal plotted over the provided time vector')

figure(1)
plot(t,eog)
grid on
xlabel('Time [sec]')
ylabel('EOG [V]')

% 9.2.b

disp('Question 9.2.b')
Mean = mean(eog);
fprintf('The mean of the signal is %d\n',Mean)

% 9.2.c

disp('Question 9.2.c')
Power = mean(eog.^2);
fprintf('The power of the signal is %d\n',Power)

% 9.2.d

disp('Question 9.2.d')
Fs = 1/(t(2)-t(1));
N = length(eog);
xdft = fft(eog);
endpoint_xdft = (N/2)+0.5;
xdft = xdft(1:endpoint_xdft);
Power_spec = abs(xdft).^2;
freq = 0:Fs/N:Fs/2;

figure(2)
plot(freq,Power_spec)
set(gca,'Xlim',[0 1])
grid on
title('Power Spectrum')
xlabel('Frequency (Hz)')
ylabel('Power')

% 9.2.e

disp('Question 9.2.e')
Max_power = max(Power_spec);
Max_freq_index = find(Power_spec == Max_power);
Max_freq = freq(Max_freq_index);

fprintf(' %d is the frequency at which the signal is most prevalent\n',Max_freq)
