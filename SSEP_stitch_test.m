%% load file
TDTPATH = 'TDTMatlabSDK';
addpath(genpath(TDTPATH));
foldername = uigetdir; 

%%
sig = RodentSSEPtoSig(foldername);
Fs = sig.SampleRate;

%% define parameters for filter 
N = 64; % filter taps 
stepsize = .2;
nUpdates = 100;

%% pre-filtering
% notch out 60Hz (& odd harmonics)
notches = [];
notchfreq = 60; 
while notchfreq < 600
notchfreq
notchw = 2;
notch60 = designfilt('bandstopiir', ...
                     'HalfPowerFrequency1', notchfreq-notchw, ...
                     'HalfPowerFrequency2', notchfreq+notchw, ...
                     'FilterOrder', 2, ...
                     'SampleRate', Fs, ...
                     'DesignMethod', 'butter');
notches = [notches, notch60];
notchfreq = notchfreq + 120;
end

% highpass filtering (baseline removal) 
hpFilt = designfilt('highpassiir', ...
                    'StopbandFrequency', 5, ...
                    'PassbandFrequency', 10, ...
                    'PassbandRipple', .5, ...
                    'StopbandAttenuation', 60, ...
                    'SampleRate', Fs, ... 
                    'DesignMethod', 'butter');

% lowpass filtering (noise removal)
lpFilt = designfilt('lowpassiir', ...
                    'StopbandFrequency', 4500, ...
                    'PassbandFrequency', 4000, ...
                    'PassbandRipple', .5, ...
                    'StopbandAttenuation', 60, ...
                    'SampleRate', Fs, ... 
                    'DesignMethod', 'cheby2');

%% signal filtering 

filts = [buildFilterObj([notches, hpFilt], lpFilt, N, stepsize, ...
                        [-1*ones(size(notches)), 0], 0, true); ...
         buildFilterObj([notches, hpFilt], lpFilt, N, stepsize, ...
                        [+1*ones(size(notches)), 0], 0, true)];
sigs = [sig; sig];

for idx = 1:length(filts)
    filti = filts(idx);
    sigi = doPreFilt(filti, sig);
    sigi = getTrainTestWrapper(sigi);
    [sigi, filti] = preTrainWtsWrapper(filti, sigi, .1*nUpdates);
    sigi = LMSonlineWrapper(filti, sigi, .1*nUpdates);
    sigi = doPostFilt(filti, sigi);
    sigs(idx) = sigi; filts(idx) = filti;
    clear sigi filti
end

%% Preview 
figure; 
plot(sigs(1).Times, sigs(1).Data_HPF); 
hold on; grid on; 
%plot(sigs(1).Times(N:end), sigs(1).Data_LMS); 
plot(sigs(1).Times(N:end), sigs(1).Data_LMS_LPF); 
xlabel('t (s)'); ylabel('Signal (V)');

%% split into pre, post 
tBeforeTrig = .02;
[tPost, dPost, ePost] = getPrePostStim(tBeforeTrig, ...
    sigs(1).Noise_Reference, sigs(1).Data_HPF, sigs(1).Data_LMS_LPF, Fs, sigs(1).Channels, N);
[tPre,  dPre,  ePre ] = getPrePostStim(tBeforeTrig, ...
    sigs(2).Noise_Reference, sigs(2).Data_HPF, sigs(2).Data_LMS_LPF, Fs, sigs(2).Channels, N);

dStitch = cell(size(dPost)); eStitch = cell(size(ePost));
for ch = 1:length(eStitch)
    dStitch{ch} = cat(3, dPre{ch}(:,:,1), dPost{ch}(:,:,2));
    eStitch{ch} = cat(3, ePre{ch}(:,:,1), ePost{ch}(:,:,2));
end

PrePostAvgAll_v2(tBeforeTrig,tPost,dStitch,eStitch,Fs,sig.Channels,10);
%PrePostAvgBatch(1,tPost,dStitch,eStitch,Fs,sig.Channels);