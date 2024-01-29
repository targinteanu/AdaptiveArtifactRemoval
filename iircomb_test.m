% Corrected example for iircomb in MATLAB

% Define sampling frequency and time vector
Fs = 1000;  % Sampling frequency in Hz
t = 0:1/Fs:1;
fnotch = 60; % Hz 

% Generate a signal with multiple harmonic components
signal = 0.1*randn(size(t));
noisy_signal = sin(2*pi*60*t) + 0.5*sin(2*pi*180*t) + 0.2*sin(2*pi*300*t) + signal;

% Specify the filter order
filter_order = Fs/(fnotch);  % Adjust as needed
filter_order = round(filter_order);

% Specify the bandwidth for notches (should be less than 1)
notch_bandwidth = 0.1;  % Adjust as needed

% Design a comb filter with multiple notches using iircomb
[numerator, denominator] = iircomb(filter_order, notch_bandwidth);

% Apply the filter to the signal using the filter function
filtered_signal = filter(numerator, denominator, noisy_signal);

% Plot the original and filtered signals
figure;
plot(t, [signal; noisy_signal; filtered_signal]); 
legend('original', 'noisy', 'filtered');
title('Filtered Signal with iircomb');

% Frequency response plot (optional)
%{
figure;
freqz(numerator, denominator);
title('Frequency Response of the iircomb Filter');
%}