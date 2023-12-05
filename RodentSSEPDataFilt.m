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
trainfrac = .1;
N = 64; % filter taps 
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

%% cleanup 
clear Dt dt_mean dt_resample dt_err tLen 
clear t2 g2 d2
clear T G dta dta_t_chan
clear g_trl t_trl t_stim chan ch chIdx 

%% pre-filtering
% highpass filtering (baseline removal) 
hpFilt = designfilt('highpassiir', ...
                    'StopbandFrequency', .5, ...
                    'PassbandFrequency', 1.5, ...
                    'PassbandRipple', .5, ...
                    'StopbandAttenuation', 60, ...
                    'SampleRate', Fs, ... 
                    'DesignMethod', 'butter');
d         = filter(hpFilt, d_unfilt);
%% lowpass filtering (noise removal)
lpFilt = designfilt('lowpassiir', ...
                    'StopbandFrequency', 2000, ...
                    'PassbandFrequency', 1000, ...
                    'PassbandRipple', .5, ...
                    'StopbandAttenuation', 60, ...
                    'SampleRate', Fs, ... 
                    'DesignMethod', 'cheby2');
%fvtool(lpFilt);

%% pretraining and testing 
splIdx = floor(trainfrac*size(t,1));
tTrainBnd = [t(1), t(splIdx)];
[w, e_train, op_train, t_train, g_train, d_train, t_test, g_test, d_test, T, G, D] = ...
    preTrainWts(t, g, d, tTrainBnd, N, uchan, .1*nUpdates);
%[e_test, op_test] = testPreTrained(w, t_test, g_test, d_test, N, uchan, .1*nUpdates);

%% online LMS 
[e_t, w_OL] = LMSonline(t, g, d, stepsize, N, uchan, w, nUpdates);

%% post-filtering
disp('LP Filtering Train Signal')
e_train_lpf = filtfilt(lpFilt, e_train);
%disp('LP Filtering Test Signal')
%e_test_lpf  = filtfilt(lpFilt, e_test);
disp('LP Filtering Online Signal')
e_t_lpf     = filtfilt(lpFilt, e_t);
disp('LP Filtering Original Signal')
d_lpf       = filtfilt(lpFilt, d);

%% demo final signal 
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

%%
%plotBeforeAfterStim(.29, g, d, e_t_lpf, Fs, uchan, N, 150, .1*nUpdates);

%%
[t_PrePost, d_PrePost, e_PrePost] = getPrePostStim(.05, g, d, e_t_lpf, Fs, uchan, N);
PrePostAvgAll(.2,t_PrePost,d_PrePost,e_PrePost,Fs,uchan,10);
PrePostAvgBatch(90,t_PrePost,d_PrePost,e_PrePost,Fs,uchan);
PrePostStats(t_PrePost,d_PrePost,e_PrePost,Fs,uchan,10);