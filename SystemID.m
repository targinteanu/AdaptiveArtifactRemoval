%% load file
TDTPATH = 'TDTMatlabSDK';
addpath(genpath(TDTPATH));
foldername = 'Rodent SSEP Data/AC5-230830-130841'; 

sig = RodentSSEPtoSig(foldername);
Fs = sig.SampleRate;

%% define parameters for filter 
N = 2048; % filter taps 
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
                    'StopbandFrequency', .5, ...
                    'PassbandFrequency', 1.5, ...
                    'PassbandRipple', .5, ...
                    'StopbandAttenuation', 60, ...
                    'SampleRate', Fs, ... 
                    'DesignMethod', 'butter');

% lowpass filtering (noise removal)
lpFilt = designfilt('lowpassiir', ...
                    'StopbandFrequency', 3000, ...
                    'PassbandFrequency', 2000, ...
                    'PassbandRipple', .5, ...
                    'StopbandAttenuation', 60, ...
                    'SampleRate', Fs, ... 
                    'DesignMethod', 'cheby2');

%% signal filtering 

filts = [buildFilterObj([notches, hpFilt], lpFilt, N, stepsize, ...
                        [-1*ones(size(notches)), 0], 0, true); ...
         buildFilterObj([notches, hpFilt], lpFilt, N, stepsize, ...
                             [+1*ones(size(notches)), 0], 0, false)];
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

%% split into pre, post 
tBeforeTrig = .06;
[tPost, dPost, ePost] = getPrePostStim(tBeforeTrig, ...
    sigs(1).Noise_Reference, sigs(1).Data_BPF, sigs(1).Data_LMS_LPF, Fs, sigs(1).Channels, N);
[tPre,  dPre,  ePre ] = getPrePostStim(tBeforeTrig, ...
    sigs(2).Noise_Reference, sigs(2).Data_BPF, sigs(2).Data_LMS_LPF, Fs, sigs(1).Channels, N);

PrePostAvgAll_v2(tBeforeTrig,tPost,dPost,ePost,Fs,sigs(1).Channels,10);
PrePostAvgAll_v2(tBeforeTrig,tPre,dPre,ePre,Fs,sigs(2).Channels,10);