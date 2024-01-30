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
% highpass filtering (baseline removal) 
hpFilt = designfilt('highpassiir', ...
                    'StopbandFrequency', .5, ...
                    'PassbandFrequency', 1.5, ...
                    'PassbandRipple', .5, ...
                    'StopbandAttenuation', 60, ...
                    'SampleRate', Fs, ... 
                    'DesignMethod', 'butter');
%d         = filtfilt(hpFilt, d_unfilt);
%% lowpass filtering (noise removal)
lpFilt = designfilt('lowpassiir', ...
                    'StopbandFrequency', 3000, ...
                    'PassbandFrequency', 2000, ...
                    'PassbandRipple', .5, ...
                    'StopbandAttenuation', 60, ...
                    'SampleRate', Fs, ... 
                    'DesignMethod', 'cheby2');
%fvtool(lpFilt);

%% filter to object 
filtObj = buildFilterObj(hpFilt, lpFilt, N, stepsize, 0, 0, true);

sig = doPreFilt(filtObj, sig);

%% pretraining and testing 
%{
[t_train, g_train, d_train, t_test, g_test, d_test] = ...
    getTrainTest(t, g, d, tTrainBnd, uchan);
[w, e_train, op_train, fig] = preTrainWts(t_train, g_train, d_train, N, uchan, .1*nUpdates);
%}
%[e_test, op_test] = testPreTrained(w, t_test, g_test, d_test, N, uchan, .1*nUpdates);

%%{
sig = getTrainTestWrapper(sig);
[sig, filtObj] = preTrainWtsWrapper(filtObj, sig, .1*nUpdates);
%}

%% online LMS 
%[e_t, w_OL] = LMSonline(t, g, d, stepsize, N, uchan, w, nUpdates);
[sig, w_end] = LMSonlineWrapper(filtObj, sig, nUpdates);

%% post-filtering
%{
disp('LP Filtering Train Signal')
e_train_lpf = filtfilt(lpFilt, e_train);
%disp('LP Filtering Test Signal')
%e_test_lpf  = filtfilt(lpFilt, e_test);
disp('LP Filtering Online Signal')
e_t_lpf     = filtfilt(lpFilt, e_t);
disp('LP Filtering Original Signal')
d_lpf       = filtfilt(lpFilt, d);
%}

sig = doPostFilt(filtObj, sig);

%% demo final signal 
%{
for idx = 1:length(uchan)
    fig = figure; 
    figure(fig); plot(t(:,idx), d(:,idx), 'k', 'LineWidth', 1); hold on;
%    figure(fig); plot(t_train(:,idx), e_train_lpf(:,idx)); plot(t_test(:,idx), e_test_lpf(:,idx));
    figure(fig); plot(t(N:end,idx), e_t_lpf(:,idx)); 
    grid on;
    xlabel('time (s)'); ylabel('filtered signal (V)');
%    legend('original', 'train', 'test', 'online');
    legend('original', 'adaptive filtered');
    title(['channel ',num2str(uchan(idx))])
    
    %xlim([1410.1, 1411.4])
    xlim([1410.351, 1410.449])
    %xlim([336.1, 337.1])
    ylim([-8e-5, 8e-5])
end

%ylim([-8e-5, 8e-5])
%xlim([336.351, 336.449])
%xlim([336.1, 337.1])
%}

%%
%plotBeforeAfterStim(.29, g, d, e_t_lpf, Fs, uchan, N, 150, .1*nUpdates);

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
for ch = 1:length(size(ERPvals,2))

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
        [n20n40,~,stat] = measureERP(tPost, d_PrePost_ch(trl,:,2), [.02 .04], [], [.01,.1], 200, false);
        tblUnfilt.n20amp(trl) = n20n40(1,1); tblUnfilt.n20lat(trl) = n20n40(2,1); % [amplitude; latency]
        tblUnfilt.n40amp(trl) = n20n40(1,2); tblUnfilt.n40lat(trl) = n20n40(2,2); % [amplitude; latency]
        tblUnfilt.mean(trl) = stat(1); tblUnfilt.SD(trl) = std(d_PrePost_ch(trl,:,2)); 
            % mean = only within selected time range
            % noise = std thru all times
        clear n20n40 stat

        %title('unfilt');
        %subplot(212); 

        % do the same for filtered data 
        [n20n40,~,stat] = measureERP(tPost, e_PrePost_ch(trl,:,2), [.02 .04], [], [.01,.1], 200, false);
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
    title(['Latency: p = ',num2str(ploc20),' (n20), ',num2str(ploc40),' (n40)']); 
    xtl = {'Unfilt n20', 'Filt n20', 'Unfilt n40', 'Filt n40'};
    xtl = arrayfun(@(i) [xtl{i},': N=',num2str(numfound(i))], 1:length(xtl), 'UniformOutput',false);
    xticklabels(xtl);
end
linkaxes(ax1,'y'); linkaxes(ax2,'y');
%}