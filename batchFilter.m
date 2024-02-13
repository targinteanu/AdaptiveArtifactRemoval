%% setup
codedir = cd;
folder = uigetdir; 
cd(folder);
files = dir('*_sigObj.mat');
cd(codedir);
sig = loadSignalObj([folder,'\',files(1).name]);
Fs = sig.SampleRate;
clear sig

%% define parameters for filter 
% to do: change this to load a saved filter from file 

N = 64; % filter taps 
stepsize = .2;
nUpdates = 100;

% notch out 60Hz (& odd harmonics)
notches = [];
notchfreq = 60; 
while notchfreq < 600
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
                    'PassbandFrequency', 2500, ...
                    'PassbandRipple', .5, ...
                    'StopbandAttenuation', 60, ...
                    'SampleRate', Fs, ... 
                    'DesignMethod', 'cheby2');

% filter to object 
filtObj = buildFilterObj([notches, hpFilt], lpFilt, N, stepsize, ...
                         [-1*ones(size(notches)), 0], 0, true);

%% loop 
savename = 'filtered_signal_obj';
if isfolder([folder,'/',savename]) | isfolder([folder,'\',savename])
    error('already saved')
else
    mkdir(folder,savename);
end
cd(folder);
savedir = [folder,'\',savename,'\'];
cd(codedir);

%%
for f = 1:length(files)
    sig = loadSignalObj([folder,'\',files(f).name]);
    if ~(sig.SampleRate == Fs)
        error('Sample rates are not equal.')
    end

    sig = doPreFilt(filtObj, sig);
    sig = getTrainTestWrapper(sig);
    [sig, filtObj] = preTrainWtsWrapper(filtObj, sig, .1*nUpdates);
    [sig, w_end] = LMSonlineWrapper(filtObj, sig, nUpdates);
    sig = doPostFilt(filtObj, sig);

    saveSignalCompact([savedir,files(f,1).name,'_filtered'],sig);
    clear sig
end