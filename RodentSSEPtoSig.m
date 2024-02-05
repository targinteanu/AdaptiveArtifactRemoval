function sig = RodentSSEPtoSig(foldername)

data = TDTbin2mat(foldername);

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

%% define parameters for training 
trainfrac = .01;

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

%% shorten - to be removed 
%{
d_unfilt = d_unfilt(1:1000000,:);
t = t(1:1000000,:);
g = g(1:1000000,:);
%}

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

% sig = extractChannel(sig, 2); % to remove 

end