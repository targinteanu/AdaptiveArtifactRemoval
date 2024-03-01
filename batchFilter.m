%% setup
codedir = cd;
folder = uigetdir; 
cd(folder);
files = dir('*_sigObj.mat');
cd(codedir);
sig = loadSignalObj([folder,filesep,files(1).name]);
Fs = sig.SampleRate;
clear sig

%% define parameters for filter - SSEP
% to do: change this to load a saved filter from file 

%{

N = 64; % filter taps 
stepsize = .2;

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

%}

%% Define parameters for filter - motor nerve
% to do: change this to load a saved filter from file 

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
stepsize = .9;

filtObj = buildFilterObj(hpFilt, lpFilt, N, stepsize, 0, 0, true, true);

%% loop 
savename = 'filtered_signal_obj';
if isfolder([folder,'/',savename]) | isfolder([folder,filesep,savename])
    error('already saved')
else
    mkdir(folder,savename);
end
cd(folder);
savedir = [folder,filesep,savename,filesep];
cd(codedir);

%%

% unsure if updates are compatible with parallel computing
%nUpdates = 100; 
nUpdates = 0;

parfor f = 1:length(files)
    disp(['Starting ',files(f).name])
    FO = filtObj;
    sig = loadSignalObj([folder,filesep,files(f).name]);
    if ~(sig.SampleRate == Fs)
        error(['in ',files(f).name,': Sample rates are not equal.'])
    end

    sig = doPreFilt(FO, sig);
    sig = getTrainTestWrapper(sig);
    [sig, FO] = preTrainWtsWrapper(FO, sig, .1*nUpdates);
    [sig, w_end] = LMSonlineWrapper(FO, sig, nUpdates);
    sig = doPostFilt(FO, sig);

    saveSignalCompact([savedir,files(f).name,'_filtered'],sig);
    % clear sig
    disp(['Completed ',files(f).name])
end