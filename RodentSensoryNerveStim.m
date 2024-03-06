filename = '220202_Experimental5Month_Sci_L_L4DRG_220202_131245_comb.rhd';

[amplifier_channels, amplifier_data, board_adc_channels,...
          board_adc_data,frequency_parameters,notes,spike_triggers,...,
          t_amplifier,t_board_adc]...
          =read_Intan_RHD2000_file(filename); 

%chansel = [2,3,6,7,8];
chansel = [2,3];
amplifier_channels = amplifier_channels(chansel);
amplifier_data = amplifier_data(chansel,:);
spike_triggers = spike_triggers(chansel);

board_adc_data = board_adc_data';
t_board_adc = t_board_adc';
amplifier_data = amplifier_data';
t_amplifier = t_amplifier';

%% timing 
fsAmp = frequency_parameters.amplifier_sample_rate;
fsADC = frequency_parameters.board_adc_sample_rate;
tdiff = mean(abs(t_amplifier-t_board_adc));
if (~(fsAmp == fsADC)) | (tdiff > 1/fsAmp) | (tdiff > 1/fsADC)
    Fs = max(fsAmp, fsADC); 
    ti = min([t_board_adc; t_amplifier]);
    tf = max([t_board_adc; t_amplifier]);
    t = ti:(1/Fs):tf; t = t';
    warning(['Resampling at ',num2str(Fs),' Hz'])
    % resample amplifier_data and board_adc_data
    amplifier_data = interp1(t_amplifier, amplifier_data, t, 'nearest','extrap');
    board_adc_data = interp1(t_board_adc, board_adc_data, t, 'nearest','extrap');
else
    t = t_amplifier; 
    Fs = fsAmp; 
end

%% processing loaded data into signal object
% "filt" is the neural recording data, which we want to change into "d_unfilt"
% "board_adc_data" is the noise reference data, which we want to change into "g"
d_unfilt = amplifier_data;

% introduce delay in signal compared to noise reference 
padT = .005; % delay dur (s); should be at least .0005
padL = floor(padT*Fs); 
dPad = d_unfilt(1,:).*ones(padL,1);
d_unfilt = [dPad; d_unfilt(1:(end-padL),:)];

d_unfilt = (d_unfilt-32768)*0.0003125; %filt in multiple of volt units

g = board_adc_data;

chans = arrayfun(@(ch) buildChannelObj(ch.custom_channel_name, 0,0,0, 'Cartesian'), amplifier_channels);

% (t and g should have the same number of rows as d_unfilt see repmat)
t = repmat(t, 1, size(d_unfilt,2));
g = repmat(g, 1, size(d_unfilt,2));

%% define parameters for filter and training 
nUpdates = 100;
splIdx = 4e5;
tTrainBnd = [t(1,1), t(splIdx,1)];

sig = buildSignalObj([], d_unfilt, t, g, Fs, chans, ...
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
N = 150; % filter taps 
stepsize = .001;
filtObj = buildFilterObj(hpFilt, lpFilt, N, stepsize, 0, 0, true, true);
sig = doPreFilt(filtObj, sig);
sig = getTrainTestWrapper(sig);
%[sig, filtObj] = preTrainWtsWrapper(filtObj, sig, .1*nUpdates);
[sig, w_end] = LMSonlineWrapper(filtObj, sig, nUpdates);
sig = doPostFilt(filtObj, sig);

%%
figure('Units','normalized', 'Position',[.1,.1,.8,.8]); 
ax(1) = subplot(211);
plot(sig.Times(:,1), sig.Data_HPF(:,1));
hold on; 
plot(sig.Times(N:end,1), sig.Data_LMS(:,1));
legend('Unfiltered', 'Filtered');
xlabel('time (s)'); ylabel('Volts'); title('signal');
grid on;
ax(2) = subplot(212); 
plot(sig.Times(:,1), sig.Noise_Reference(:,1));
xlabel('time (s)'); ylabel('current (amps?)'); title('noise reference');
grid on;
linkaxes(ax,'x');