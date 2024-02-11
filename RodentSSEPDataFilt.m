%% load file
TDTPATH = 'TDTMatlabSDK';
addpath(genpath(TDTPATH));
foldername = 'Rodent SSEP Data/AC5-230830-130841'; 

sig = RodentSSEPtoSig(foldername);
Fs = sig.SampleRate;

%% define parameters for filter 
N = 64; % filter taps 
stepsize = .2;
nUpdates = 100;

%% pre-filtering
% notch out 60Hz (& odd harmonics)
notches = [];
%%{
notchfreq = 60; 
while notchfreq < 600
notchfreq
%notchw = notchfreq/40;
notchw = 2;
%%{
notch60 = designfilt('bandstopiir', ...
                     'HalfPowerFrequency1', notchfreq-notchw, ...
                     'HalfPowerFrequency2', notchfreq+notchw, ...
                     'FilterOrder', 2, ...
                     'SampleRate', Fs, ...
                     'DesignMethod', 'butter');
%}
%[notchN,notchD] = iirnotch(notchfreq/(Fs/2), notchw/(Fs/2));
%fvtool(notch60);
notches = [notches, notch60];
notchfreq = notchfreq + 120;
end
%}

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

%% filter to object 
filtObj = buildFilterObj([notches, hpFilt], lpFilt, N, stepsize, ...
                         [-1*ones(size(notches)), 0], 0, true);

%% pre-filtering 
sig = doPreFilt(filtObj, sig);

%% pretraining 
sig = getTrainTestWrapper(sig);
[sig, filtObj] = preTrainWtsWrapper(filtObj, sig, .1*nUpdates);

%% online LMS 
[sig, w_end] = LMSonlineWrapper(filtObj, sig, nUpdates);

%% post-filtering
sig = doPostFilt(filtObj, sig);

%% Preview 
figure; 
plot(sig.Times, sig.Data_HPF); 
hold on; grid on; 
%plot(sig.Times(N:end), sig.Data_LMS); 
plot(sig.Times(N:end), sig.Data_LMS_LPF); 
xlabel('t (s)'); ylabel('Signal (V)');

%%
tBeforeTrig = .06;
[t_PrePost, d_PrePost, e_PrePost] = getPrePostStim(tBeforeTrig, ...
    sig.Noise_Reference, sig.Data_BPF, sig.Data_LMS_LPF, Fs, sig.Channels, N);
PrePostAvgAll_v2(tBeforeTrig,t_PrePost,d_PrePost,e_PrePost,Fs,sig.Channels,10);

%{
PrePostAvgAll(tBeforeTrig,t_PrePost,d_PrePost,e_PrePost,Fs,sig.Channels,10);
PrePostAvgBatch(90,t_PrePost,d_PrePost,e_PrePost,Fs,sig.Channels);
PrePostStats(t_PrePost,d_PrePost,e_PrePost,Fs,sig.Channels,[1.5,1000],10);
%}

%% SNR (signal : noise ratio)
tPost = t_PrePost(2,:);
ERPvals = cell(2, length(e_PrePost)); % rows = unfilt vs filt; cols = channel

% calculate values for all channels 
for ch = 1:size(ERPvals,2)

    % get filtered and unfiltered waveforms at this channel for each trial
    e_PrePost_ch = e_PrePost{ch};
    d_PrePost_ch = d_PrePost{ch};

    % get number of trials 
    nTrl = size(d_PrePost_ch,1);

    % setup tables for calculated values 
    tblFilt = table('Size',[nTrl,6],...
                    'VariableTypes',["double","double","double","double","double","double"], ...
                    'VariableNames',{'n20amp','n20lat','n40amp','n40lat', 'mean','SD'});
    tblUnfilt = tblFilt;

    % calculate values for each trial 
    for trl = 1:nTrl

        %figure; subplot(211); 

        % populate row trl of table with n20 & n40 amps and latencies and stats (mean & SD)
        [n20n40,~,stat] = measureERP(tPost, d_PrePost_ch(trl,:,2), [.015 .007], [], [.004,.1], 200, false);
        tblUnfilt.n20amp(trl) = n20n40(1,1); tblUnfilt.n20lat(trl) = n20n40(2,1); % [amplitude; latency]
        tblUnfilt.n40amp(trl) = n20n40(1,2); tblUnfilt.n40lat(trl) = n20n40(2,2); % [amplitude; latency]
        tblUnfilt.mean(trl) = stat(1); tblUnfilt.SD(trl) = std(d_PrePost_ch(trl,:,2)); 
            % mean = only within selected time range
            % noise = std thru all times
        clear n20n40 stat

        %title('unfilt');
        %subplot(212); 

        % do the same for filtered data 
        [n20n40,~,stat] = measureERP(tPost, e_PrePost_ch(trl,:,2), [.015 .007], [], [.004,.1], 200, false);
        tblFilt.n20amp(trl) = n20n40(1,1); tblFilt.n20lat(trl) = n20n40(2,1); % [amplitude; latency]
        tblFilt.n40amp(trl) = n20n40(1,2); tblFilt.n40lat(trl) = n20n40(2,2); % [amplitude; latency]
        tblFilt.mean(trl) = stat(1); tblFilt.SD(trl) = std(e_PrePost_ch(trl,:,2)); 
        clear n20n40 stat

        %title('filt');
    end

    % store tables in cell array 
    ERPvals{1, ch} = tblUnfilt; ERPvals{2, ch} = tblFilt;
    clear tblFilt tblUnfilt e_PrePost_ch d_PrePost_ch;
end

%%
figure('Units','normalized', 'Position',[.1,.1,.8,.8])
for ch = 1:size(ERPvals,2)
    % get tables for this channel 
    tblUnfilt = ERPvals{1,ch}; tblFilt = ERPvals{2,ch}; 
    % "Signal" = n20 amplitude - mean for [unfilt, filt]
    SNR_S = [tblUnfilt.n20amp, tblFilt.n20amp] - [tblUnfilt.mean, tblFilt.mean];
    % "Noise" = SD for [unfilt, filt]
    SNR_N = [tblUnfilt.SD, tblFilt.SD];
    % SNR = S / N
    SNR = SNR_S./SNR_N; 
    % hypothesis test whether unfilt has greater SNR than filt 
    [~,pSNR] = ttest(SNR(:,1), SNR(:,2), 'Tail', 'left');
    % loc = n20, n40 latency for unfilt, filt 
    loc = [tblUnfilt.n20lat, tblFilt.n20lat, tblUnfilt.n40lat, tblFilt.n40lat];
    numfound = sum(~isnan(loc));
    % hypothesis test whether unfilt and filt have different locs 
    [~,ploc20] = ttest(loc(:,1), loc(:,2), 'Tail', 'both');
    [~,ploc40] = ttest(loc(:,3), loc(:,4), 'Tail', 'both');

    % box plot results with p values displayed  
    ax1(ch) = subplot(size(ERPvals,2),2,2*(ch-1)+1); boxplot(SNR);
    grid on;
    title(['SNR: p = ',num2str(pSNR)]); ylabel(['Channel ',sig.Channels(ch).labels]);
    xticklabels({'Unfilt', 'Filt'});
    ax2(ch) = subplot(size(ERPvals,2),2,2*ch); boxplot(loc);
    grid on;
    title(['Latency: p = ',num2str(ploc20),' (n20), ',num2str(ploc40),' (n7)']); 
    xtl = {'Unfilt n20', 'Filt n20', 'Unfilt n7', 'Filt n7'};
    xtl = arrayfun(@(i) [xtl{i},': N=',num2str(numfound(i))], 1:length(xtl), 'UniformOutput',false);
    xticklabels(xtl);
end
linkaxes(ax1,'y'); linkaxes(ax2,'y');
%}