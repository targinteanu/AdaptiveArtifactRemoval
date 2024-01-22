%% load file
load('rat1/amplifier_data/Naive/rat1_N_CMAP_1_230731_113111'); % load the "filt" 2D array 
%load('D:\filtering research proj\rat1\rat1\amplifier_data\Naive\rat1_N_CMAP_1_230731_113111.mat');
load('rat1/trigger_data/Naive/rat1_N_CMAP_1_230731_113111_trigger'); % load the "board_adc_data" 1D array 
%load('D:\filtering research proj\rat1\rat1\trigger_data\Naive\rat1_N_CMAP_1_230731_113111_trigger.mat');

%% define timing
Fs = 20000; % samples per second 
dt = 1/Fs; % time step (seconds) 
samples=1:length(filt);
samples = samples-1;
t = samples/Fs;
t = t';

%% processing loaded data into signal object
% "filt" is the neural recording data, which we want to change into "d_unfilt"
% "board_adc_data" is the noise reference data, which we want to change into "g"
d_unfilt = filt;

% introduce delay in signal compared to noise reference 
padT = .005; % delay dur (s); should be at least .0005
padL = floor(padT*Fs); 
dPad = d_unfilt(:,1).*ones(1,padL);
d_unfilt = [dPad, d_unfilt(:,1:(end-padL))];

d_unfilt = (d_unfilt-32768)*0.0003125; %filt in multiple of volt units
d_unfilt = d_unfilt';

%g = [board_adc_data;board_adc_data];
g = board_adc_data;
g = g';

chA = buildChannelObj('Chan1', 0,0,0, 'Cartesian');
chB = buildChannelObj('Chan2', 0,0,0, 'Cartesian');

% (t and g should have the same number of rows as d_unfilt see repmat)
t = repmat(t, 1, size(d_unfilt,2));
g = repmat(g, 1, size(d_unfilt,2));
%% define parameters for filter and training 
trainfrac = .1;
nUpdates = 100;
splIdx = floor(trainfrac*size(t,1));
tTrainBnd = [t(1), t(splIdx)];

sig = buildSignalObj([], d_unfilt, t, g, Fs, [chA; chB], ...
                     tTrainBnd, tTrainBnd, 2);

%% filtering
% highpass filtering (baseline removal) 
hpFilt = designfilt('highpassiir', ...
                    'StopbandFrequency', .5, ...
                    'PassbandFrequency', 1.5, ...
                    'PassbandRipple', .5, ...
                    'StopbandAttenuation', 60, ...
                    'SampleRate', Fs, ... 
                    'DesignMethod', 'butter');
% lowpass filtering (noise removal)
lpFilt = designfilt('lowpassiir', ...
                    'StopbandFrequency', 1500, ...
                    'PassbandFrequency', 1000, ...
                    'PassbandRipple', .5, ...
                    'StopbandAttenuation', 60, ...
                    'SampleRate', Fs, ... 
                    'DesignMethod', 'cheby2');
N = 512; % filter taps 
stepsize = .01;
filtObj = buildFilterObj(hpFilt, lpFilt, N, stepsize, true);
sig = doHPFilt(filtObj, sig);
%sig = getTrainTestWrapper(sig);
%[sig, filtObj] = preTrainWtsWrapper(filtObj, sig, .1*nUpdates);
[sig, w_end] = LMSonlineWrapper(filtObj, sig, nUpdates);
sig = doLPFilt(filtObj, sig);

%%
plotwindow1 = 1:1100000;
plotwindow2 = plotwindow1 + N;
figure('Units','normalized', 'Position',[.1,.1,.8,.8]); 
ax(1) = subplot(211);
%plot(sig.Times(plotwindow1), sig.Data_BPF(plotwindow1,2));
plot(sig.Times(plotwindow1), sig.Data_HPF(plotwindow1,2));
hold on; 
%plot(sig.Times(plotwindow2), sig.Data_LMS_LPF(plotwindow1,2));
plot(sig.Times(plotwindow2), sig.Data_LMS(plotwindow1,2));
legend('Unfiltered', 'Filtered');
xlabel('time (s)'); ylabel('Volts'); title('signal');
grid on;
ax(2) = subplot(212); 
plot(sig.Times(plotwindow1), sig.Noise_Reference(plotwindow1,2));
xlabel('time (s)'); ylabel('current (amps?)'); title('noise reference');
grid on;
linkaxes(ax,'x');