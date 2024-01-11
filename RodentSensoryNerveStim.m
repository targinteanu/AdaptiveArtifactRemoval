%% load file
%load('rat1/amplifier_data/Naive/rat1_N_CMAP_1_230731_113111'); % load the "filt" 2D array 
load('D:\filtering research proj\rat1\rat1\amplifier_data\Naive\rat1_N_CMAP_1_230731_113111.mat');
%load('rat1/trigger_data/Naive/rat1_N_CMAP_1_230731_113111_trigger'); % load the "board_adc_data" 1D array 
load('D:\filtering research proj\rat1\rat1\trigger_data\Naive\rat1_N_CMAP_1_230731_113111_trigger.mat');

%% define timing--no response yet from Kiara's team?
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
d_unfilt = (d_unfilt-32768)*0.0003125; %filt in multiple of volt units
d_unfilt = d_unfilt';

g = [board_adc_data;board_adc_data];
g = g';

chA = buildChannelObj('insert_name_here', 0,0,0, 'Cartesian');
chB = buildChannelObj('insert_name_here', 0,0,0, 'Cartesian');

% (t and g should have the same number of rows as d_unfilt see repmat)
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
                    'StopbandFrequency', 2000, ...
                    'PassbandFrequency', 1500, ...
                    'PassbandRipple', .5, ...
                    'StopbandAttenuation', 60, ...
                    'SampleRate', Fs, ... 
                    'DesignMethod', 'cheby2');
N = 128; % filter taps 
stepsize = .1;
filtObj = buildFilterObj(hpFilt, lpFilt, N, stepsize, true);
sig = doHPFilt(filtObj, sig);
[sig, w_end] = LMSonlineWrapper(filtObj, sig, nUpdates);


plotwindow = 300000:600000;
figure; 
plot(sig.Times(plotwindow), sig.Data_HPF(plotwindow,2));
hold on; 
plot(sig.Times(plotwindow), sig.Data_LMS(plotwindow,2));
legend('Unfiltered', 'Filtered');
grid on;