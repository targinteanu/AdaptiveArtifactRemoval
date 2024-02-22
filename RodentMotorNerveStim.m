%% load files
%{
load('rat1/amplifier_data/Naive/rat1_N_CMAP_1_230731_113111'); % load the "filt" 2D array 
%load('D:\filtering research proj\rat1\rat1\amplifier_data\Naive\rat1_N_CMAP_1_230731_113111.mat');
load('rat1/trigger_data/Naive/rat1_N_CMAP_1_230731_113111_trigger'); % load the "board_adc_data" 1D array 
%load('D:\filtering research proj\rat1\rat1\trigger_data\Naive\rat1_N_CMAP_1_230731_113111_trigger.mat');
%}
foldername = 'rat1/amplifier_data/Naive';
files = dir(foldername);
files = files(~[files.isdir]);
d0 = []; g0 = [];
for f = 1:length(files)
    files(f).name
    load([foldername,'/',files(f).name]);
    d0 = [d0, filt];
end
foldername = 'rat1/trigger_data/Naive';
files = dir(foldername);
files = files(~[files.isdir]);
for f = 1:length(files)
    files(f).name
    load([foldername,'/',files(f).name]);
    g0 = [g0, board_adc_data];
end

clear filt board_adc_data

%%{
figure; 
subplot(211); plot(d0');
subplot(212); plot(g0);
%}

%% truncate signal 
%d0 = d0(:,1:3450000);
%g0 = g0(:,1:3450000);
% d0 = d0(:,1:500000);
% g0 = g0(:,1:500000);

%% define timing
Fs = 20000; % samples per second 
dt = 1/Fs; % time step (seconds) 
samples=1:length(d0);
samples = samples-1;
t = samples/Fs;
t = t';

%% fix noise reference 
% each amplitude: 10 pulses 
ampvals = {20:2:120, ...   % 2 uA => 20uA – 120uA
           125:5:300, ...  % 5 uA => 125uA – 300uA
           310:10:500, ... % 10 uA => 310uA – 500uA
           550:50:750};    % 50 uA => 550uA - 750uA
ampvalSz = arrayfun(@(i) length(ampvals{i}), 1:length(ampvals))';
ampvals = cell2mat(ampvals);

g = g0; % reset 
g = g./max(g);
g = g.*(g > .1);
% find actual recorded number of pulses 
[trigval, trigloc] = findpeaks(g, 'MinPeakProminence', .1*max(g));

% split pulses into groups of 10 
nInGrp = 10;
ngrp = length(trigloc)/nInGrp;
pulseFirstLast = zeros(nInGrp,ngrp);
pulseFirstLast(:) = trigloc(:);
pulseFirstLast = pulseFirstLast';
pulseFirstLast = pulseFirstLast(:,[1,end]);
%{
pulseFirstLast = zeros(ngrp, 2);
[k,C] = kmeans(trigloc', ngrp);
[Cord,ord] = sort(C); ku = unique(k); kuord = ku(ord);
for idx = 1:ngrp
    kIdx = kuord(idx); Cidx = Cord(idx);
    trigidx = trigloc(k==kIdx);
    pulseFirstLast(idx,1) = min(trigidx); pulseFirstLast(idx,2) = max(trigidx);
end
%}
pulsesGap0 = [pulseFirstLast(2:end,1), pulseFirstLast(1:(end-1),2)];
pulsesBnd = mean(pulsesGap0, 2)'; 

if sum(ampvalSz) == ngrp
    % There are the same number of expected and recorded pulses. Proceed
    % assigning amplitude from beginning. 
    startFromEnd = false; 
else
    % There is a mismatch between expected vs recorded number of pulses.
    % Decide how to assign amplitudes. 

pulsesGap = diff(pulsesGap0, [], 2);
[~,ord] = sort(pulsesGap);
ord = ord(1:3); ord = sort(ord);
groupFirstLast = [[0;ord]+1, [ord;ngrp]];
%groupSz = diff(groupFirstLast, [], 2);
groupSz = zeros(size(groupFirstLast,1),1);

for idx = 1:size(groupFirstLast,1)
    groupIdx = groupFirstLast(idx,:);
    pulsesIdx = pulseFirstLast(groupIdx(1):groupIdx(2), :);
    trigidx = (trigloc >= pulsesIdx(1,1)) & (trigloc <= pulsesIdx(end,end));
    groupSz(idx) = sum(trigidx)/10;
end

mismatchSz = ampvalSz - groupSz; 
if (mismatchSz(1) <= 0) & (mismatchSz(end) <= 0)
    % error is in the middle? 
    error('Unable to resolve discrepancy between expected and recorded stimuli.')
elseif mismatchSz(1) <= 0
    % end is cut off but beginning is not. 
    startFromEnd = false;
elseif mismatchSz(end) <= 0
    % beginning is cut off but end is not. 
    startFromEnd = true;
else
    % mismatch at both beginning and end
    if abs(mismatchSz(1)) <= abs(mismatchSz(end))
        startFromEnd = false;
        ampvals = ampvals(abs(mismatchSz(1)):end);
    else
        startFromEnd = true;
        ampvals = ampvals(1:(end-abs(mismatchSz(end))));
    end
end

end

if ~startFromEnd
    bnd1 = 1; idx = 1;
    pulsesBnd = [pulsesBnd, length(g)];
    for bnd2 = pulsesBnd
        g(bnd1:bnd2) = g(bnd1:bnd2) * ampvals(idx);
        idx = idx + 1; bnd1 = bnd2+1;
    end
else
    bnd2 = length(g); idx = 1; 
    ampvals = fliplr(ampvals);
    pulsesBnd = [1, pulsesBnd];
    for bnd1 = fliplr(pulsesBnd)
        g(bnd1:bnd2) = g(bnd1:bnd2) * ampvals(idx);
        idx = idx + 1; bnd2 = bnd1-1;
    end
end

figure; plot(g)

%% processing loaded data into signal object
% "filt" is the neural recording data, which we want to change into "d_unfilt"
% "board_adc_data" is the noise reference data, which we want to change into "g"
d_unfilt = d0;

% introduce delay in signal compared to noise reference 
padT = .005; % delay dur (s); should be at least .0005
padL = floor(padT*Fs); 
dPad = d_unfilt(:,1).*ones(1,padL);
d_unfilt = [dPad, d_unfilt(:,1:(end-padL))];

d_unfilt = (d_unfilt-32768)*0.0003125; %filt in multiple of volt units
d_unfilt = d_unfilt';

%g = [board_adc_data;board_adc_data];
%g = board_adc_data;
g = g';

chA = buildChannelObj('Chan1', 0,0,0, 'Cartesian');
chB = buildChannelObj('Chan2', 0,0,0, 'Cartesian');

% (t and g should have the same number of rows as d_unfilt see repmat)
t = repmat(t, 1, size(d_unfilt,2));
g = repmat(g, 1, size(d_unfilt,2));
%% define parameters for filter and training 
%trainfrac = .02;
nUpdates = 100;
%splIdx = floor(trainfrac*size(t,1));
splIdx = 14e4;
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
                    'StopbandFrequency', 1500, ...
                    'PassbandFrequency', 1000, ...
                    'PassbandRipple', .5, ...
                    'StopbandAttenuation', 60, ...
                    'SampleRate', Fs, ... 
                    'DesignMethod', 'cheby2');
N = 150; % filter taps 
stepsize = .9;
filtObj = buildFilterObj(hpFilt, lpFilt, N, stepsize, 0, 0, true, true);
sig = doPreFilt(filtObj, sig);
sig = getTrainTestWrapper(sig);
[sig, filtObj] = preTrainWtsWrapper(filtObj, sig, .1*nUpdates);
[sig, w_end] = LMSonlineWrapper(filtObj, sig, nUpdates);
sig = doPostFilt(filtObj, sig);

%%
%plotwindow1 = 1:1100000;
plotwindow2 = 3450000 - ((1e6):(-1):0);
plotwindow1 = plotwindow2 - N;
%plotwindow2 = plotwindow1 + N;
figure('Units','normalized', 'Position',[.1,.1,.8,.8]); 
ax(1) = subplot(211);
%plot(sig.Times(plotwindow1), sig.Data_BPF(plotwindow1,2));
plot(sig.Times(plotwindow1), sig.Data_HPF(plotwindow1,2));
hold on; 
%plot(sig.Times(plotwindow2), sig.Data_LMS_LPF(plotwindow1,2));
plot(sig.Times(plotwindow2), sig.Data_LMS(plotwindow1,2));
legend('Unfiltered', 'Filtered');
xlabel('time (s)'); ylabel('Volts'); title('signal');
grid on;
ax(2) = subplot(212); 
plot(sig.Times(plotwindow1), sig.Noise_Reference(plotwindow1,2));
xlabel('time (s)'); ylabel('current (amps?)'); title('noise reference');
grid on;
linkaxes(ax,'x');

%%
tBeforeTrig = .02;
[t_PrePost, d_PrePost, e_PrePost] = getPrePostStim(tBeforeTrig, ...
    sig.Noise_Reference, sig.Data_BPF, sig.Data_LMS_LPF, Fs, sig.Channels, N);
PrePostAvgAll_v2(tBeforeTrig,t_PrePost,d_PrePost,e_PrePost,Fs,sig.Channels,10);

%%{
%PrePostAvgAll(tBeforeTrig,t_PrePost,d_PrePost,e_PrePost,Fs,sig.Channels,10);
PrePostAvgBatch(10,t_PrePost,d_PrePost,e_PrePost,Fs,sig.Channels);
%PrePostStats(t_PrePost,d_PrePost,e_PrePost,Fs,sig.Channels,[1.5,1000],10);
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
                    'VariableNames',{'p10amp','p10lat','n14amp','n14lat', 'mean','SD'});
    tblUnfilt = tblFilt;

    % calculate values for each trial 
    for trl = 1:nTrl

        %figure; subplot(211); 

        % populate row trl of table with n20 & n40 amps and latencies and stats (mean & SD)
        [n14,p10,stat] = measureERP(tPost, d_PrePost_ch(trl,:,2), .014, .01, [.008,.25], 10, false);
        p10n14 = [p10, n14];
        tblUnfilt.p10amp(trl) = p10n14(1,1); tblUnfilt.p10lat(trl) = p10n14(2,1); % [amplitude; latency]
        tblUnfilt.n14amp(trl) = p10n14(1,2); tblUnfilt.n14lat(trl) = p10n14(2,2); % [amplitude; latency]
        tblUnfilt.mean(trl) = stat(1); tblUnfilt.SD(trl) = std(d_PrePost_ch(trl,:,2)); 
            % mean = only within selected time range
            % noise = std thru all times
        clear p10n14 stat p10 n14

        %hold on; plot(tPost, d_PrePost_ch(trl,:,2));
        %title('unfilt');
        %subplot(212); 

        % do the same for filtered data 
        [n14,p10,stat] = measureERP(tPost, e_PrePost_ch(trl,:,2), .014, .01, [.008,.25], 10, false);
        p10n14 = [p10, n14];
        tblFilt.p10amp(trl) = p10n14(1,1); tblFilt.p10lat(trl) = p10n14(2,1); % [amplitude; latency]
        tblFilt.n14amp(trl) = p10n14(1,2); tblFilt.n14lat(trl) = p10n14(2,2); % [amplitude; latency]
        tblFilt.mean(trl) = stat(1); tblFilt.SD(trl) = std(e_PrePost_ch(trl,:,2)); 
        clear p10n14 stat p10 n14

        %hold on; plot(tPost, e_PrePost_ch(trl,:,2));
        %title('filt');
    end

    % store tables in cell array 
    ERPvals{1, ch} = tblUnfilt; ERPvals{2, ch} = tblFilt;
    clear tblFilt tblUnfilt e_PrePost_ch d_PrePost_ch;
end

%%
figure('Units','normalized', 'Position',[.05,.1,.9,.8])
for ch = 1:size(ERPvals,2)
    % get tables for this channel 
    tblUnfilt = ERPvals{1,ch}; tblFilt = ERPvals{2,ch}; 
    % "Signal" = n14 - p10 amplitude for [unfilt, filt]
    SNR_S = [tblUnfilt.n14amp, tblFilt.n14amp] - [tblUnfilt.p10amp, tblFilt.p10amp];
    % "Noise" = SD for [unfilt, filt]
    SNR_N = [tblUnfilt.SD, tblFilt.SD];
    % SNR = S / N
    SNR = SNR_S./SNR_N; 
    % hypothesis test whether unfilt has greater SNR than filt 
    [~,pSNR] = ttest(SNR(:,1), SNR(:,2), 'Tail', 'left');
    % amp = p10, n14 amplitude for unfilt, filt 
    amp = [tblUnfilt.p10amp, tblFilt.p10amp, tblUnfilt.n14amp, tblFilt.n14amp];
    % loc = p10, n14 latency for unfilt, filt 
    loc = [tblUnfilt.p10lat, tblFilt.p10lat, tblUnfilt.n14lat, tblFilt.n14lat];
    numfound = sum(~isnan(loc));
    % hypothesis test whether unfilt and filt have different locs, amps 
    [~,ploc10] = ttest(loc(:,1), loc(:,2), 'Tail', 'both');
    [~,ploc14] = ttest(loc(:,3), loc(:,4), 'Tail', 'both');
    [~,pamp10] = ttest(amp(:,1), amp(:,2), 'Tail', 'both');
    [~,pamp14] = ttest(amp(:,3), amp(:,4), 'Tail', 'both');

    % box plot results with p values displayed  
    ax1(ch) = subplot(size(ERPvals,2),3,3*(ch-1)+1); boxplot(SNR);
    grid on;
    title(['SNR: p = ',num2str(pSNR)]); ylabel(['Channel ',sig.Channels(ch).labels]);
    xticklabels({'Unfilt', 'Filt'});
    ax2(ch) = subplot(size(ERPvals,2),3,3*ch-1); boxplot(amp);
    grid on; 
    title(['Amplitude: p = ',num2str(pamp10),' (p10), ',num2str(pamp14),' (n14)']); 
    xtl = {'Unfilt p10', 'Filt p10', 'Unfilt n14', 'Filt n14'};
    xtl = arrayfun(@(i) [xtl{i},': N=',num2str(numfound(i))], 1:length(xtl), 'UniformOutput',false);
    xticklabels(xtl);
    ax3(ch) = subplot(size(ERPvals,2),3,3*ch); boxplot(loc);
    grid on;
    title(['Latency: p = ',num2str(ploc10),' (p10), ',num2str(ploc14),' (n14)']); 
    %xtl = {'Unfilt p10', 'Filt p10', 'Unfilt n14', 'Filt n14'};
    %xtl = arrayfun(@(i) [xtl{i},': N=',num2str(numfound(i))], 1:length(xtl), 'UniformOutput',false);
    xticklabels(xtl);
end
linkaxes(ax1,'y'); linkaxes(ax2,'y');
%}

%% Sys ID 
splIdx = round(length(t)/2);
tTrainBnd = [t(splIdx-5e5), t(splIdx+5e5)];
sig.Train_Time_Bounds = tTrainBnd;
sig = getTrainTestWrapper(sig);
[TTtrain,TTtest,TT] = getTrainTestTimetable(sig, filtObj); 
sysLin = n4sid(TTtrain,4, 'OutputName',["y1","y2"],'InputName','u'); 
sysTF = tfest(TTtrain,2, 'OutputName',["y1","y2"],'InputName','u');
figure; compare(TTtrain, sysLin);
figure; compare(TT, sysLin);
sysHW0 = nlhw(TTtrain,sysLin); 
figure; compare(TTtrain, sysHW0); 
sysGB = idnlgrey(TTtrain,sysLin);
% seems like a good idea but takes too long 
%{
sigmoidIn = idSigmoidNetwork(1); waveletOut = idWaveletNetwork(1);
sysHW = nlhw(TTtrain, sysLin, sigmoidIn, waveletOut);
figure; compare(TTtrain, sysHW); 
%}
% try idSigmoidNetwork input, idWaveletNetwork output 
% network of 1 faster?
%figure; compare(TT, sysHW);