%% load file
TDTPATH = 'TDTMatlabSDK';
addpath(genpath(TDTPATH));
data = TDTbin2mat('Rodent SSEP Data/AC5-230830-130841');

%% define parameters of stim 
% monophasic, fixed
%  ^         __________
%  |        |          |
% AmpA      |          |
%  |        |          |
%  v   _____|          |_____
%            <--DurA-->
PeriodA = 0.8;  % ms - between pulses?
CountA  = 1;    % pulse count? 
AmpA    = 1500; % uA
DurA    = 0.2;  % ms
DelayA  = 0;    % ms
ChanA   = 1; 

%% organize data from file
dta = data.snips.SSEP.data; 
t_stim = data.snips.SSEP.ts;
fs = data.snips.SSEP.fs;
chan = data.snips.SSEP.chan; uchan = unique(chan);

dta_t_chan = cell(2, length(uchan));
for idx = 1:length(uchan)
    ch = uchan(idx);
    chIdx = chan == ch;
    dta_t_chan{1, idx} =  dta(chIdx,:);
    dta_t_chan{2, idx} = t_stim(chIdx);
end

dta = zeros([size(dta_t_chan{1,1}), length(uchan)]);
t_stim = zeros([length(dta_t_chan{2,1}), length(uchan)]);
for idx = 1:length(uchan)
    dta(:,:,idx)  = dta_t_chan{1, idx};
    t_stim(:,idx) = dta_t_chan{2, idx};
end

t_trl = (1:size(dta, 2))/fs - .3; % ~ -.3 to +1 s
g_trl = (AmpA/1000)*((t_trl >= 0)&(t_trl < DurA/1000)); % noise reference, mA
G = repmat(g_trl, size(dta,1), 1, size(dta,3));

T = zeros(size(dta));
for idx = 1:length(uchan)
    T(:,:,idx) = t_stim(:,idx) + t_trl;
end

%% define parameters for filter and training 
trainfrac = .01;
N = 1024; % filter taps 
stepsize = .2;
nUpdates = 100;

%% "linearize" trial blocks 
%uchan = uchan(1); % comment out to get all chans
t        = zeros(size(T,1)  *size(T,2),   length(uchan));
g        = zeros(size(G,1)  *size(G,2),   length(uchan));
d_unfilt = zeros(size(dta,1)*size(dta,2), length(uchan));
for idx = 1:length(uchan)
    Tidx = T(:,:,idx)'; Gidx = G(:,:,idx)'; Didx = dta(:,:,idx)';

    % ensure correct order of timepoints 
    [Tidx, ord] = sort(Tidx(:));

    t(ord,idx)        = Tidx(:);
    g(ord,idx)        = Gidx(:);
    d_unfilt(ord,idx) = Didx(:);
end

% detect and fix inconsistencies in sampling 
% t must be in columns!
Dt = diff(t);
dt_mean = mean(Dt(:));
dt_resample  = 1/fs;
dt_err  = std(Dt(:));
tLen = t(end,:) - t(1,:);
if dt_err > .01*dt_mean
    warning(['Inconsistent time steps; resampling at ',num2str(1/dt_resample),' Hz']);
    t_from_start = 0:dt_resample:max(tLen);

    t2 = zeros(length(t_from_start), length(uchan));
    g2 = zeros(length(t_from_start), length(uchan));
    d2 = zeros(length(t_from_start), length(uchan));

    for idx = 1:length(uchan)
        if tLen(idx) < max(tLen)
            warning(['Channel ',num2str(uchan(idx)),' has shorter duration and may be end-padded']);
        end
        t_ch = t_from_start + t(1,idx);

        t2(:,idx) = t_ch;
        g2(:,idx) = interp1(t(:,idx),        g(:,idx), t_ch, 'nearest','extrap');
        d2(:,idx) = interp1(t(:,idx), d_unfilt(:,idx), t_ch, 'nearest','extrap');
    end
    t        = t2; 
    g        = g2; 
    d_unfilt = d2; 

    dt_mean = dt_resample;
end

Fs = 1/dt_mean; % Hz

splIdx = floor(trainfrac*size(t,1));
tTrainBnd = [t(1), t(splIdx)];

%% notch out 60Hz (& odd harmonics?)
% +- ~15
% to do
notch60 = designfilt('bandstopiir', ...
                     'PassbandFrequency1', 45, ...
                     'StopbandFrequency1', 55, ...
                     'StopbandFrequency2', 65, ...
                     'PassbandFrequency2', 75, ...
                     'PassbandRipple1', .5, ...
                     'PassbandRipple2', .5, ...
                     'StopbandAttenuation', 20, ...
                     'SampleRate', Fs, ...
                     'DesignMethod', 'cheby2');
%fvtool(notch60);
d_unfilt_2 = filtfilt(notch60, d_unfilt);
d_unfilt = d_unfilt_2;

%% shorten - to be removed 
d_unfilt = d_unfilt(1:2000000,:);
t = t(1:2000000,:);
g = g(1:2000000,:);

%% cleanup 
%{
clear Dt dt_mean dt_resample dt_err tLen 
clear t2 g2 d2
clear T G dta dta_t_chan
clear g_trl t_trl t_stim chan ch chIdx 
%}

%% signal to object 
chA = buildChannelObj('A', 1, 1,0,'Cartesian'); % right somatosensory
chB = buildChannelObj('B',-1, 1,0,'Cartesian'); % left somatosensory
chC = buildChannelObj('C', 1,-1,0,'Cartesian'); % ?
chD = buildChannelObj('D',-1,-1,0,'Cartesian'); % ?

sig = buildSignalObj([], d_unfilt, t, g, Fs, [chA; chB; chC; chD], ...
                     tTrainBnd, tTrainBnd, 2);

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
                    'StopbandFrequency', 6000, ...
                    'PassbandFrequency', 5000, ...
                    'PassbandRipple', .5, ...
                    'StopbandAttenuation', 60, ...
                    'SampleRate', Fs, ... 
                    'DesignMethod', 'cheby2');
%fvtool(lpFilt);

%% filter to object 
filtObj = buildFilterObj(hpFilt, lpFilt, N, stepsize, true, true);

sig = doHPFilt(filtObj, sig);

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

sig = doLPFilt(filtObj, sig);

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
tBeforeTrig = .2;
[t_PrePost, d_PrePost, e_PrePost] = getPrePostStim(tBeforeTrig, ...
    sig.Noise_Reference, sig.Data_BPF, sig.Data_LMS_LPF, Fs, sig.Channels, N);
PrePostAvgAll_v2(tBeforeTrig,t_PrePost,d_PrePost,e_PrePost,Fs,sig.Channels,10);

%{
PrePostAvgAll(tBeforeTrig,t_PrePost,d_PrePost,e_PrePost,Fs,sig.Channels,10);
PrePostAvgBatch(90,t_PrePost,d_PrePost,e_PrePost,Fs,sig.Channels);
PrePostStats(t_PrePost,d_PrePost,e_PrePost,Fs,sig.Channels,[1.5,1000],10);

%% SNR 
tPost = t_PrePost(2,:);
n20s = cell(size(e_PrePost));
for ch = 1:length(n20s)
    e_PrePost_ch = e_PrePost{ch};
    d_PrePost_ch = d_PrePost{ch};
    nTrl = size(d_PrePost_ch,1);
    n20 = zeros(3,nTrl,2);
    for trl = 1:nTrl
        n20([1,2],trl,1) = measureERP(tPost, d_PrePost_ch(trl,:,2), .02, [], [.01,.05]);
        n20(3,trl,1) = std(d_PrePost_ch(trl,:,2));
        n20([1,2],trl,2) = measureERP(tPost, e_PrePost_ch(trl,:,2), .02, [], [.01,.05]);
        n20(3,trl,2) = std(e_PrePost_ch(trl,:,2));
    end
    n20s{ch} = n20;
    clear n20 e_PrePost_ch d_PrePost_ch;
end

%%
figure('Units','normalized', 'Position',[.1,.1,.8,.8])
for ch = 1:length(n20s)
    n20 = n20s{ch};
    SNR = n20(1,:,:)./n20(3,:,:); SNR = squeeze(SNR);
    [~,pSNR] = ttest(SNR(:,1), SNR(:,2), 'Tail', 'left');
    loc = n20(2,:,:); loc = squeeze(loc);
    [~,ploc] = ttest(loc(:,1), loc(:,2), 'Tail', 'both');
    % titles and axes 
    ax1(ch) = subplot(length(n20s),2,2*(ch-1)+1); boxplot(SNR);
    grid on;
    title(['SNR: p = ',num2str(pSNR)]); ylabel(['Channel ',sig.Channels(ch).labels]);
    xticklabels({'Unfilt', 'Filt'});
    ax2(ch) = subplot(length(n20s),2,2*ch); boxplot(loc);
    grid on;
    title(['Latency: p = ',num2str(ploc)]); 
    xticklabels({'Unfilt', 'Filt'});
end
linkaxes(ax1,'y'); linkaxes(ax2,'y');
%}