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

%{
figure; 
subplot(211); plot(d0);
subplot(212); plot(g0);
%}

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
trainfrac = .1;
nUpdates = 100;
splIdx = floor(trainfrac*size(t,1));
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
N = 512; % filter taps 
stepsize = .000001;
filtObj = buildFilterObj(hpFilt, lpFilt, N, stepsize, 0, 0, true);
sig = doPreFilt(filtObj, sig);
sig = getTrainTestWrapper(sig);
[sig, filtObj] = preTrainWtsWrapper(filtObj, sig, .1*nUpdates);
[sig, w_end] = LMSonlineWrapper(filtObj, sig, nUpdates);
sig = doPostFilt(filtObj, sig);

%%
plotwindow1 = 1:1100000;
plotwindow2 = plotwindow1 + N;
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
tBeforeTrig = .03;
[t_PrePost, d_PrePost, e_PrePost] = getPrePostStim(tBeforeTrig, ...
    sig.Noise_Reference, sig.Data_BPF, sig.Data_LMS_LPF, Fs, sig.Channels, N);
PrePostAvgAll_v2(tBeforeTrig,t_PrePost,d_PrePost,e_PrePost,Fs,sig.Channels,10);

%%{
%PrePostAvgAll(tBeforeTrig,t_PrePost,d_PrePost,e_PrePost,Fs,sig.Channels,10);
PrePostAvgBatch(10,t_PrePost,d_PrePost,e_PrePost,Fs,sig.Channels);
%PrePostStats(t_PrePost,d_PrePost,e_PrePost,Fs,sig.Channels,[1.5,1000],10);
%}

%% SNR 
tPost = t_PrePost(2,:);
p10s = cell(size(e_PrePost));
n14s = cell(size(p10s));
for ch = 1:length(p10s)
    e_PrePost_ch = e_PrePost{ch};
    d_PrePost_ch = d_PrePost{ch};
    nTrl = size(d_PrePost_ch,1);
    p10 = zeros(3,nTrl,2); n14 = p10;
    for trl = 1:nTrl
        [p10([1,2],trl,1), n14([1,2],trl,1)] = ...
            measureERP(tPost, d_PrePost_ch(trl,:,2), .01, .014, [.007,.25]);
        p10(3,trl,1) = std(d_PrePost_ch(trl,:,2));
        [p10([1,2],trl,2), n14([1,2],trl,2)] = ...
            measureERP(tPost, e_PrePost_ch(trl,:,2), .01, .014, [.007,.25]);
        p10(3,trl,2) = std(e_PrePost_ch(trl,:,2));
    end
    p10s{ch} = p10; n14s{ch} = n14;
    clear p10 n14 e_PrePost_ch d_PrePost_ch;
end

%%
figure('Units','normalized', 'Position',[.1,.1,.8,.8])
for ch = 1:length(p10s)
    p10 = p10s{ch}; n14 = n14s{ch};
    SNR = (n14(1,:,:) - p10(1,:,:))./p10(3,:,:); SNR = squeeze(SNR);
    [~,pSNR] = ttest(SNR(:,1), SNR(:,2), 'Tail', 'left');
    loc = n14(2,:,:); loc = squeeze(loc);
    [~,ploc] = ttest(loc(:,1), loc(:,2), 'Tail', 'both');
    % titles and axes 
    ax1(ch) = subplot(length(p10s),2,2*(ch-1)+1); boxplot(SNR);
    grid on;
    title(['SNR: p = ',num2str(pSNR)]); ylabel(['Channel ',sig.Channels(ch).labels]);
    xticklabels({'Unfilt', 'Filt'});
    ax2(ch) = subplot(length(p10s),2,2*ch); boxplot(loc);
    grid on;
    title(['Latency: p = ',num2str(ploc)]); 
    xticklabels({'Unfilt', 'Filt'});
end
linkaxes(ax1,'y'); linkaxes(ax2,'y');
%}

%% Sys ID 
t = sig.Times(N:end,1); u = sig.Noise_Reference(N:end,1); 
y1 = sig.Data_LMS(:,1); y2 = sig.Data_LMS(:,2); 
TT = timetable(seconds(t),u,y1,y2);
sys = ssest(TT,'InputName','u','OutputName',["y1","y2"]);